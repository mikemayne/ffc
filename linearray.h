#include "ffc.h"
#include "frd/frd.h"

// LineArray defines an array of homogeneous elements and models
// their mechanical properties and rigging.
// LineArray provides the method createSoundSources()
// to make the data structure the ffc needs
template<typename FloatType>
struct LineArray
{
    LineArray(int numElements) :
    polar (nullptr),
    elementHeight (0.3),
    elementDepth  (0.4),
    validSplayAngles_ ({FloatType(0)}),
    inverted (false),
    arrayAimAngle {0.0}
    { 
        for (int i=0; i<numElements; ++i)
            elements.push_back({validSplayAngles_[0]});
    }

    ffc::Polar<FloatType> const* polar;
    FloatType elementHeight;
    FloatType elementDepth;
    FloatType arrayAimAngle;
    std::pair<FloatType, FloatType> arrayXYinM;

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
        Element (FloatType splayAngle) :
        splayAngle(splayAngle)
        { }

        FloatType splayAngle; // this enclosures splay angle
        FloatType cumulativeAngle; // the difference between vertical and this enclosure.
        std::pair<FloatType, FloatType> topRiggingXY;
        std::pair<FloatType, FloatType> bottomRiggingXY;
        const FloatType height;
        const FloatType width;
    };
    std::vector<Element> elements;

    const std::vector<FloatType>& validSplayAngles() const {
        return validSplayAngles_;
    }

    // void setup_linearray(std::pair<FloatType, FloatType> xy)
    // {
    //     auto angleAccumulator = arrayAimAngle; 
    //     auto hangPoint = xy;
    //     auto isFlown = inverted ? static_cast<FloatType>(-1) : static_cast<FloatType>(1);

    //     for (auto& element : elements) 
    //     {
    //         angleAccumulator += element.splayAngle;
    //         element.source.angle = angleAccumulator;

    //         auto h = element.height;
    //         auto dx = isFlown * h * sin (degreesToRadians (angleAccumulator));
    //         auto dy = isFlown * h * cos (degreesToRadians (angleAccumulator));

    //         element.topRiggingXY = hangPoint;
    //         element.acousticCentre = {hangPoint.first + dx/2., hangPoint.second + dy/2.};

    //         hangPoint.first += dx;
    //         hangPoint.second += dy;
    //         element.bottomRiggingXY = hangPoint;
    //     }
    // }

    // 
    void recalc(std::pair<FloatType, FloatType> originXY) {
        const auto toRadians = M_PI / 180.0;
        FloatType angleAccumulator = arrayAimAngle;
        std::pair<FloatType, FloatType> xyAccumulator = {FloatType(0), FloatType(0)};
        for (auto &e : elements) {
            e.topRiggingXY = xyAccumulator;
            angleAccumulator += e.splayAngle;
            e.cumulativeAngle = angleAccumulator;
            xyAccumulator.first -= elementHeight * std::sin(angleAccumulator * toRadians);
            xyAccumulator.second -= elementHeight * std::cos(angleAccumulator * toRadians);
            std::cout << "xyAccumulator " << xyAccumulator.first << ", " << xyAccumulator.second << " angle " << angleAccumulator << "deg" << std::endl;
            e.bottomRiggingXY = xyAccumulator;
        }
    }

    // creates the data structure needed by the calculator
    // assumes the acoustic centre is in the middle of the baffle
    std::vector<ffc::SoundSource<FloatType>> createSoundSources()
    {
        std::vector<ffc::SoundSource<FloatType>> soundSources;
        auto angleAccumulator = FloatType(0);
        for (const auto& e : elements)
        {
            FloatType x = e.topRiggingXY.first  + e.bottomRiggingXY.first  / FloatType(2);
            FloatType y = e.topRiggingXY.second + e.bottomRiggingXY.second / FloatType(2);

            angleAccumulator += e.splayAngle;

            // angle is going to be an issue, need to accumulate them through the array
            soundSources.push_back(ffc::SoundSource<FloatType>::create({x, y}, polar, angleAccumulator));
        }
        return soundSources;
    }
private:
    std::vector<FloatType> validSplayAngles_;
    bool inverted;
};
