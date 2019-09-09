#pragma once

#include <vector>
#include <math.h>
#include <utility>      // std::pair
#include <complex>		// std::complex
#include <algorithm>
#include <functional>
#include <thread>

#include "frd/frd.h"

namespace ffc 
{
template <typename FloatType>
FloatType to_degrees(FloatType radians)
{
    return radians * 180 / M_PI;
}

template <typename FloatType>
struct Constants
{
    FloatType frequencyHz;
    FloatType m_per_pixel;
    FloatType speed_of_sound_mps;
    int width_px;
    int height_px;
    std::pair<FloatType, FloatType> XYM; // origin in M (bad member name!)
};

// Polar holds magnitude and phase in RADIANS for a single frequency.  
// Its used in the calculators critical path
template<typename FloatType>
using Polar = std::vector<std::pair<FloatType, FloatType>>;

template<typename FloatType>
FloatType constrainAngle (FloatType x) {
    x = fmod(x + 180,360);
    if (x < 0)
        x += 360;
    return x - 180;
}

template<typename FloatType>
Polar<FloatType> make_polar(frd::PolarData<FloatType> const& polarData, FloatType frequencyHz)
{
    auto toRadians = M_PI / 180.0;

    Polar<FloatType> polar;
    polar.reserve(polarData.size());
    for (auto frdList : polarData) {
        auto value = frd::findFreq(frdList, frequencyHz);
        if (value != frdList.end()) {
            auto valueToInsert = std::make_pair(value->dBSPL, value->phaseDeg * toRadians);
            polar.push_back(valueToInsert);
        } else {
            throw std::runtime_error("Couldn't find frequency " + std::to_string(frequencyHz));
        }
    }
    return polar;
}

template<typename FloatType>
Polar<FloatType> make_mirrored_polar(frd::PolarData<FloatType> const& polarData, FloatType frequencyHz)
{
    Polar<FloatType> polar = make_polar<FloatType>(polarData, frequencyHz);
    polar.reserve (2*polar.size());
    Polar<FloatType> reverse {polar.rbegin()+1, polar.rend()-1};
    polar.insert(polar.end(), reverse.begin(), reverse.end());
    return polar;
}

template<typename FloatType>
std::pair<FloatType, FloatType> interpolated(Polar<FloatType> const& p, FloatType rads)
{
    rads += M_PI;
    auto rads_per_element = 2 * (FloatType)M_PI / (FloatType)p.size();
    auto step_below = static_cast<typename Polar<FloatType>::size_type>(floor(rads / rads_per_element));

    auto index_below = step_below % p.size();
    auto index_above = (index_below + 1) % p.size();

    auto t = (rads - step_below * rads_per_element) / rads_per_element; // 0 < t < 1

    auto spl_phase_below = p[index_below];
    auto spl_phase_above = p[index_above];

    auto spl   = spl_phase_below.first  + t * (spl_phase_above.first  - spl_phase_below.first);
    auto phase = spl_phase_below.second + t * (spl_phase_above.second - spl_phase_below.second);

    return std::make_pair(spl, phase);
}

template<typename FloatType>
std::tuple<FloatType, FloatType, FloatType> calc_intermediate(const std::pair<FloatType, FloatType> acousticCentreXYM, const std::pair<FloatType, FloatType> XYM, const Constants<FloatType>& constants)
{
    auto dx_dy = std::make_pair((acousticCentreXYM.first - XYM.first), (acousticCentreXYM.second - XYM.second));
    auto distance_m = hypot(dx_dy.first, dx_dy.second);
    distance_m = distance_m < 1 ? 1 : distance_m; // minimum distance 1m

    auto t_secs = distance_m / constants.speed_of_sound_mps;
    auto wt = 2 * M_PI * constants.frequencyHz * t_secs;
    auto angle_rads = atan2(dx_dy.second, dx_dy.first);
    auto attenuation_dB = 20 * log10(distance_m);

    return std::make_tuple(attenuation_dB, wt, angle_rads);
}

template<typename FloatType>
FloatType convertSplToMag (FloatType spl) {
    return pow (10.0, spl / 20.0); 
}

template<typename FloatType>
FloatType convertMagtoSpl(std::complex<FloatType> mag) {
    return 20 * log10 (abs(mag));
}

template<typename FloatType>
std::vector<FloatType> convert_to_spl(std::vector<std::complex<FloatType>> const& results)
{
    std::vector<FloatType> spl(results.size());
    std::transform(results.begin(), results.end(), spl.begin(), [](std::complex<FloatType> r) -> FloatType { return convertMagtoSpl(r); });
    return spl;
}

// SoundSource is the acoustic model of the loudspeaker used by the calculator
// Its an immutable POD type,
template<typename FloatType>
struct SoundSource
{
    static SoundSource<FloatType> create  ( std::pair<FloatType, FloatType> acousticCentre,
                                            Polar<FloatType> const* polar,
                                            FloatType angle) {
        return SoundSource<FloatType>(acousticCentre, polar, angle, 1.f);
    }

    static SoundSource<FloatType> createReversePolarity (  std::pair<FloatType, FloatType> acousticCentre,
                                                Polar<FloatType> const* polar,
                                                FloatType angle) {
        return SoundSource<FloatType>(acousticCentre, polar, angle, -1.f);
    }

    const std::pair<FloatType, FloatType> acousticCentre;
    Polar<FloatType> const* polar;
    const FloatType angle;
    const FloatType polarity;

private:
    explicit SoundSource(std::pair<FloatType, FloatType> acousticCentre,
                         Polar<FloatType> const* polar,
                         FloatType angle,
                         FloatType polarity) :
                         acousticCentre (acousticCentre),
                         polar (polar),
                         angle (angle),
                         polarity (polarity) {}
};

// LineArray defines an array of homogeneous elements and models
// their mechanical properties and rigging.
// LineArray provides the method createSoundSources()
// to make the data structure the ffc needs
template<typename FloatType>
struct LineArray
{
    LineArray() :
    polar (nullptr),
    elementHeight (0.3),
    elementDepth  (0.4),
    validSplayAngles_ ({FloatType(0)}),
    inverted (false)
    { }

    Polar<FloatType> const* polar;
    FloatType elementHeight;
    FloatType elementDepth;

    bool setValidSplayAngles(std::vector<FloatType> angles) {
        if (angles.size() == 0) return false;
        validSplayAngles_ = angles;
        return true;
    }

    // non-copyable (use boost::noncopyable?)
    LineArray<FloatType>& operator=(const LineArray<FloatType>&) = delete;
    LineArray<FloatType>(const LineArray<FloatType>&) = delete;

    // Element models the state of rigging hardware etc
    struct Element
    {
        Element (LineArray<FloatType> const& lineArray) :
                splayAngle(lineArray.validSplayAngles()[0])
        { }

        FloatType splayAngle;
        std::pair<FloatType, FloatType> topRiggingXY;
        std::pair<FloatType, FloatType> bottomRiggingXY;
        struct Geometry {
            std::pair<FloatType, FloatType> rotationPoint;
            std::pair<FloatType, FloatType> xy;
            std::pair<FloatType, FloatType> hw;
            FloatType angleSum;
        };
    };
    std::vector<Element> elements;

    const std::vector<FloatType>& validSplayAngles() const {
        return validSplayAngles_;
    }

    std::vector<typename Element::Geometry> createGeometry() const {
        std::vector<typename Element::Geometry> gs;
        FloatType angleAccumulator = 0;
        for (const auto &e : elements) {
            typename Element::Geometry g;
            g.xy = {e.topRiggingXY.first - elementDepth,
                    e.topRiggingXY.second - (inverted ? elementHeight : 0.f)};
            g.hw = {elementDepth, elementHeight};
            g.rotationPoint = inverted ? e.topRiggingXY : std::make_pair(e.topRiggingXY.first, e.topRiggingXY.second + elementHeight);
            if (inverted) {
                angleAccumulator += e.splayAngle;
                g.angleSum = angleAccumulator;
            } else {
                g.angleSum = angleAccumulator;
                angleAccumulator += e.splayAngle;
            }
            gs.push_back(g);
        }
        return gs;
    }

    // creates the data structure needed by the calculator
    std::vector<SoundSource<FloatType>> createSoundSources()
    {
        std::vector<SoundSource<FloatType>> soundSources;
        auto angleAccumulator = FloatType(0);
        for (const auto& e : elements)
        {
            FloatType x = e.topRiggingXY.first  + e.bottomRiggingXY.first  / FloatType(2);
            FloatType y = e.topRiggingXY.second + e.bottomRiggingXY.second / FloatType(2);

            angleAccumulator += e.splayAngle;

            // angle is going to be an issue, need to accumulate them through the array
            soundSources.push_back(SoundSource<FloatType>::create({x, y}, polar, angleAccumulator));
        }
        return soundSources;
    }
private:
    std::vector<FloatType> validSplayAngles_;
    bool inverted;
};

template <typename FloatType>                                                   
std::pair<FloatType, FloatType> operator+(const std::pair<FloatType, FloatType> & l, const std::pair<FloatType, FloatType> & r) {   
    return {l.first+r.first,l.second+r.second};                                    
}

template<typename FloatType>
FloatType degreesToRadians(FloatType degrees) { 
    return degrees * M_PI / (FloatType)180;
}

template<typename FloatType>
struct FreeFieldCalculator
{
    FreeFieldCalculator (int numThreads, std::vector<SoundSource<FloatType>> const& sources, Constants<FloatType> const& constants)
    : stop (false)
    , numThreads(numThreads)
    , sources (sources)
    , constants (constants)
    { }

    ~FreeFieldCalculator () { }

    typedef typename std::vector<std::complex<FloatType>> Result;

    // int x_y_to_index(int x, int y, int width) {
    //     return y * width + x;
    // }

    // index_to_xy returns x,y in metres from the origin
    std::pair<FloatType, FloatType> index_to_xym(long index, Constants<FloatType> const& c) {
        auto y_px =  + index / c.width_px;
        auto x_px = index % c.width_px;
        return std::make_pair(c.XYM.first + x_px * c.m_per_pixel, c.XYM.second + y_px * c.m_per_pixel);
    }

    // calculate will block until calculation complete
    // calculate splits a vector up into equal non-overlapping segments
    // and then calculates a result in a number of threads
    void calculate(Result& result)
    {
        std::vector <std::thread> threadPool(numThreads);
        auto range = result.size() / numThreads;
        
        typename Result::iterator begin, end = result.begin();
        for (int i=0; i<numThreads; ++i) {
            begin = end;
            end = std::next(begin, range);
            threadPool.push_back(std::thread(std::bind(&FreeFieldCalculator::calculateInThread, this, result.begin(), begin, end)));
        }

        for (size_t i=0; i<threadPool.size(); ++i)
        {
            try {
                if (threadPool[i].joinable())
                    threadPool[i].join();
            } catch (std::exception const& e) {
//                std::cerr << "exception " << e.what() << " occurred joining thread " << i << std::endl;
            }
            
        }
    }
private:
    // calculateInThread is intended to be used by multiple threads to
    // calculate results in a single vector.  If its used with an
    // overlapping range of a single vector, bad things may happen.
    // concurrent modification of a std::vector is thread safe as long
    // as each index is only modified by one thread
    void calculateInThread(typename Result::iterator start, typename Result::iterator startOfRange, typename Result::iterator endOfRange)
    {
        for (const auto& s : sources) {
            for (auto pt = startOfRange; pt != endOfRange && !stop; ++pt) {
                auto xym = index_to_xym(std::distance(start, pt), constants);
                auto [attenuation_dB, wt, angle_rads] = calc_intermediate(s.acousticCentre, xym, constants);
                auto spl_phase = interpolated(*s.polar, angle_rads + s.angle);
                auto mag = convertSplToMag(spl_phase.first - attenuation_dB);
                auto phase = spl_phase.second + wt;
                *pt = *pt + std::polar(mag, phase) * s.polarity;
            }
        }
    }

    std::atomic<bool> stop;
    int numThreads;
    std::vector<SoundSource<FloatType>> const& sources;
    Constants<FloatType> const& constants;
};

}
