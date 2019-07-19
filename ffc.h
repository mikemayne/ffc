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
    std::pair<FloatType, FloatType> XYM;
};

// Polar holds magnitude and phase in RADIANS for a single frequency.  
// Its used in the calculators critical path
template<typename FloatType>
using Polar = std::vector<std::pair<FloatType, FloatType>>;

template<typename FloatType>
struct Source
{
    Source(Polar<FloatType> const& polar)
    : acousticCentre({0.f, 0.f}), polar(polar), polarity(1), angle(0)
    { }

    std::pair<FloatType, FloatType> acousticCentre;

    void reversePolarity()
    {
        polarity = -1;
    }
    Polar<FloatType> const& polar;
    FloatType polarity;
    FloatType angle;
};

template<typename FloatType>
FloatType constrainAngle (FloatType x) {
    x = fmod(x + 180,360);
    if (x < 0)
        x += 360;
    return x - 180;
}

template<typename FloatType>
Polar<FloatType> make_polar(frd::PolarData<FloatType> polarData, FloatType frequencyHz)
{
    auto toRadians = M_PI / 180.0;
    Polar<FloatType> polar(polarData.size());
    for (auto frdList : polarData) {
        if (auto value = findFreq(frdList, frequencyHz) != frdList.end()) {
            polar.push_back(std::make_pair(value->dBSPL, value->phaseDeg * toRadians));
        } else {
            throw std::runtime_error("Couldn't find frequency " + std::to_string(frequencyHz));
        }
    }
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

template<typename FloatType>
struct Element
{
    Element(Polar<FloatType> const& polar, FloatType height, FloatType depth, std::vector<FloatType> validSplayAngles)
    : source({polar}), splayAngle(FloatType(0)), splayInverted (false),
    height(height), depth(depth), validSplayAngles(validSplayAngles)
    { }

    FloatType splayAngle;
    std::pair<FloatType, FloatType> topRiggingXY;
    std::pair<FloatType, FloatType> bottomRiggingXY;
    bool splayInverted;

    const FloatType height;
    const FloatType depth;
    const std::vector<FloatType> validSplayAngles;
    Source<FloatType> source; 
};

template <typename FloatType>                                                   
std::pair<FloatType, FloatType> operator+(const std::pair<FloatType, FloatType> & l, const std::pair<FloatType, FloatType> & r) {   
    return {l.first+r.first,l.second+r.second};                                    
}

template<typename FloatType>
FloatType degreesToRadians(FloatType degrees) { 
    return degrees * M_PI / (FloatType)180;
}
// setup_linearray          set the angles and acoustic centres for a vector of elements
// xy                   the point to hang the array.
// arrayAimAngle      the initial splay angle of the array
// lineArray            the vector of elements representing the array
template<typename FloatType>
void setup_linearray(std::pair<FloatType, FloatType> xy, FloatType arrayAimAngle, std::vector<Element<FloatType>>& lineArray, bool flown)
{
    auto angleAccumulator = arrayAimAngle; 
    auto hangPoint = xy;
    auto isFlown = flown ? static_cast<FloatType>(-1) : static_cast<FloatType>(1);

    for (auto& element : lineArray) 
	{
		angleAccumulator += element.splayAngle;
		element.source.angle = angleAccumulator;

		auto h = element.height;
        auto dx = isFlown * h * sin (degreesToRadians (angleAccumulator));
        auto dy = isFlown * h * cos (degreesToRadians (angleAccumulator));

		element.topRiggingXY = hangPoint;
		element.source.acousticCentre = {hangPoint.first + dx/2., hangPoint.second + dy/2.};

    	hangPoint.first += dx;
    	hangPoint.second += dy;
    	element.bottomRiggingXY = hangPoint;
    }
}

template<typename FloatType>
struct FreeFieldCalculator
{
    FreeFieldCalculator (int numThreads, std::vector<Source<FloatType>> const& sources, Constants<FloatType> const& constants) 
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
    std::pair<FloatType, FloatType> index_to_xym(int index, Constants<FloatType> const& c) {
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
                std::cerr << "exception " << e.what() << " occurred joining thread " << i << std::endl;
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
                auto spl_phase = interpolated(s.polar, angle_rads + s.angle);
                auto mag = convertSplToMag(spl_phase.first - attenuation_dB);
                auto phase = spl_phase.second + wt;
                *pt = *pt + std::polar(mag, phase) * s.polarity;
            }
        }
    }

    std::atomic<bool> stop;
    int numThreads;
    std::vector<Source<FloatType>> const& sources;
    Constants<FloatType> const& constants;
};

}
