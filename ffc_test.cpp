#include <gtest/gtest.h>
#include "ffc.h"
#include "frd/frd.h"
#include <future>
#include "linearray.h"

namespace {
    TEST(ModelTest, MakePolar) {
        std::stringstream ss;
        for (int i=0; i<32; ++i)
            ss << i << ".0 10.0 0.0\n";
        auto frd0 = frd::parse<float>(ss);

        frd::PolarData<float> testPolarData = {frd0, frd0};
        EXPECT_EQ(testPolarData.size(), 2);
        
        auto polar = ffc::make_polar(testPolarData, 2.f);
        EXPECT_EQ(polar.size(), 2);

        auto value = ffc::interpolated(polar, 0.f);
        EXPECT_NEAR(value.first, 10.f, 0.01f);
        EXPECT_NEAR(value.second, 0.f, 0.01f);

        auto mirror = ffc::make_mirrored_polar<float>(testPolarData, 2.f);
        EXPECT_EQ(2*(polar.size()-1), mirror.size());

    }

    TEST(ModelTest, Polar) {
        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({spl, -3.f}); // minus180deg
        p.push_back({spl, -1.f}); // minus90deg
        p.push_back({spl, 1.f});  // onaxis
        p.push_back({spl, 3.f});  // plus90deg

        auto onaxis = ffc::interpolated(p, 0.f);
        EXPECT_NEAR(onaxis.second, 1.f, 0.01f);
        
        auto plus90deg = ffc::interpolated(p, (float)M_PI / 2.f);
        EXPECT_NEAR(plus90deg.second, 3.f, 0.01f); 

        auto minus180deg = ffc::interpolated(p, (float)M_PI);
        EXPECT_NEAR(minus180deg.second, -3.f, 0.01f);        
        
        auto minus90deg = ffc::interpolated(p, -(float)M_PI / 2.f);
        EXPECT_NEAR(minus90deg.second, -1.f, 0.01f);
    }

    TEST(ModelTest, Attenuation) {
        ffc::Constants<float> c;
        c.m_per_pixel = 1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 10;
        c.height_px = 10;
    
        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({98.f, 0.f});

        ffc::SoundSource<float> s = ffc::SoundSource<float>::create({0.f, 0.f}, &p, 0.f);
        std::vector<ffc::SoundSource<float>> sources = {s};
        ffc::FreeFieldCalculator<float> calc(1, sources, c);
        ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px); 
        calc.calculate(result);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[1]), 0.1f);
        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[10]), 0.1f);
        EXPECT_NEAR(spl-6.02, ffc::convertMagtoSpl(result[2]), 0.1f);
        EXPECT_NEAR(spl-6.02, ffc::convertMagtoSpl(result[20]), 0.1f);
        EXPECT_NEAR(spl-12.04, ffc::convertMagtoSpl(result[4]), 0.1f);
        EXPECT_NEAR(spl-12.04, ffc::convertMagtoSpl(result[40]), 0.1f);
    }

    TEST(ModelTest, Intermediate) {
        ffc::Constants<float> c;
        c.m_per_pixel = 1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;

        auto [att0, wt0, rads0] = calc_intermediate({0.f, 0.f}, {1.f, 0.f}, c);
        auto [att1, wt1, rads1] = calc_intermediate({1.f, 0.f}, {0.f, 0.f}, c);

        EXPECT_NEAR(att0, att1, 0.001f);    
        EXPECT_NEAR(wt0, wt1, 0.001f);        
        EXPECT_NEAR(rads0, (float)M_PI, 0.001f);        
        EXPECT_NEAR(rads1, 0.f, 0.001f);

        auto [att2, wt2, radsAbove] = calc_intermediate({0.f, 0.f}, {0.f, -1.f}, c);
        EXPECT_NEAR(radsAbove, (float)M_PI / 2.f, 0.001f);        

        auto [att3, wt3, radsBelow] = calc_intermediate({0.f, 0.f}, {0.f, 1.f}, c);
        EXPECT_NEAR(radsBelow, -(float)M_PI / 2.f, 0.001f);        
    }

    TEST(ModelTest, Sum) {
        ffc::Constants<float> c;
        c.m_per_pixel = 1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 6;
        c.height_px = 1;

        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({98.f, 0.f});
        ffc::SoundSource<float> s1 = ffc::SoundSource<float>::create({0.f, 0.f}, &p, 0.f);

        {
        std::vector<ffc::SoundSource<float>> sources = {s1};
        ffc::FreeFieldCalculator<float> calc(1, sources, c);
        ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
        calc.calculate(result);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[1]), 0.1f);
        EXPECT_NEAR(spl-6.02, ffc::convertMagtoSpl(result[2]), 0.1f);
        }

        ffc::SoundSource<float> s2 = ffc::SoundSource<float>::create({4.f, 0.f}, &p, 0.f);
        std::vector<ffc::SoundSource<float>> sources = {s1, s2};
        ffc::FreeFieldCalculator<float> calc(1, sources, c);
        ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
        calc.calculate(result);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[2]), 0.1f);
    }

    TEST(ModelTest, Cancel) {
        ffc::Constants<float> c;
        c.m_per_pixel = 1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 6;
        c.height_px = 1;

        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({98.f, 0.f});
        ffc::SoundSource<float> s1 = ffc::SoundSource<float>::create({0.f, 0.f}, &p, 0.f);

        std::vector<ffc::SoundSource<float>> sources1 = {s1};
        ffc::FreeFieldCalculator<float> calc1(2, sources1, c);
        ffc::FreeFieldCalculator<float>::Result result1(c.height_px * c.width_px);
        calc1.calculate(result1);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result1[1]), 0.1f);
        EXPECT_NEAR(spl-6.02, ffc::convertMagtoSpl(result1[2]), 0.1f);

        ffc::SoundSource<float> s2 = ffc::SoundSource<float>::createReversePolarity({4.f, 0.f}, &p, 0.f);

        std::vector<ffc::SoundSource<float>> sources2 = {s2};
        ffc::FreeFieldCalculator<float> calc2(3, sources2, c);
        ffc::FreeFieldCalculator<float>::Result result2(c.height_px * c.width_px);
        calc2.calculate(result2);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result2[3]), 0.1f);
        EXPECT_NEAR(spl-6.02, ffc::convertMagtoSpl(result2[2]), 0.1f);

        EXPECT_TRUE(result1[2] == -1.f * result2[2]);

        EXPECT_TRUE(std::isinf(ffc::convertMagtoSpl(result1[2] + result2[2])));
        EXPECT_TRUE(ffc::convertMagtoSpl(result1[2] + result2[2]) < 0.f);

        EXPECT_EQ(result1.size(), result2.size());
        EXPECT_EQ(result1.size(), c.width_px);

        std::vector<ffc::SoundSource<float>> sources3 = {s1, s2};
        ffc::FreeFieldCalculator<float> calc3(4, sources3, c);
        ffc::FreeFieldCalculator<float>::Result result3(c.height_px * c.width_px);
        calc3.calculate(result2);

        EXPECT_TRUE(std::isinf(ffc::convertMagtoSpl(result3[2])));
        EXPECT_TRUE(ffc::convertMagtoSpl(result3[2]) < 0.f);
    }

    TEST(ModelTest, Aim) {
        ffc::Constants<float> c;
        c.m_per_pixel = 1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 3;
        c.height_px = 3;

        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({spl, 0.f});
        for (int i=0; i<35; ++i) p.push_back({1.f, 0.f});

        auto s = ffc::SoundSource<float>::create({1.f, 1.f}, &p, 0.f);
        {
        std::vector<ffc::SoundSource<float>> sources = {s};
        ffc::FreeFieldCalculator<float> calc(4, sources, c);
        ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
        calc.calculate(result);

        EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[5]), 0.1f);

        EXPECT_TRUE(result[0] == result[6]);
        EXPECT_TRUE(result[1] == result[7]);
        EXPECT_NEAR(abs(result[2]), abs(result[8]), 1.f); // this is magnitude not dB so 1.f is negligible
        EXPECT_TRUE(result[3] != result[5]);
        }
        auto s_pi_2 = ffc::SoundSource<float>::create({1.f, 1.f}, &p, M_PI / 2.f);
        {
            std::vector<ffc::SoundSource<float>> sources = {s_pi_2};
            ffc::FreeFieldCalculator<float> calc(4, sources, c);
            ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
            calc.calculate(result);

            EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[1]), 0.1f);
        }
        auto s_pi = ffc::SoundSource<float>::create({1.f, 1.f}, &p, M_PI);
        {
            std::vector<ffc::SoundSource<float>> sources = {s_pi};
            ffc::FreeFieldCalculator<float> calc(4, sources, c);
            ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
            calc.calculate(result);
            EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[3]), 0.1f);

        }
        auto s_3_pi_2 = ffc::SoundSource<float>::create({1.f, 1.f}, &p, 3.f*M_PI / 2.f);
        {
            std::vector<ffc::SoundSource<float>> sources = {s_3_pi_2};
            ffc::FreeFieldCalculator<float> calc(4, sources, c);
            ffc::FreeFieldCalculator<float>::Result result(c.height_px * c.width_px);
            calc.calculate(result);
            EXPECT_NEAR(spl, ffc::convertMagtoSpl(result[7]), 0.1f);

        }
    }
    
    TEST(ModelTest, m_per_pixel_ne_1) {
        ffc::Constants<float> c;
        c.m_per_pixel = 0.1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 60;
        c.height_px = 1;
        c.XYM = {0,0};

        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({98.f, 0.f});
        
        ffc::SoundSource<float> s = ffc::SoundSource<float>::create({0.f, 0.f}, &p, 0.f);

        std::vector<ffc::SoundSource<float>> sources = {s};
        ffc::FreeFieldCalculator<float> calc(10, sources, c);
        ffc::FreeFieldCalculator<float>::Result r(c.height_px * c.width_px);
        calc.calculate(r);

        EXPECT_NEAR(spl,       ffc::convertMagtoSpl(r[10]), 0.1f);
        EXPECT_NEAR(spl-6.02,  ffc::convertMagtoSpl(r[20]), 0.1f);
        EXPECT_NEAR(spl-12.04, ffc::convertMagtoSpl(r[40]), 0.1f);
    }

    TEST(ModelTest, convert_to_spl) {
        ffc::Constants<float> c;
        c.m_per_pixel = 0.1;
        c.speed_of_sound_mps = 343;
        c.frequencyHz = 1000;
        c.width_px = 60;
        c.height_px = 1;

        ffc::Polar<float> p;
        auto spl = 98.f;
        p.push_back({98.f, 0.f});
        ffc::SoundSource<float> s = ffc::SoundSource<float>::create({0.f, 0.f}, &p, 0.f);

        std::vector<ffc::SoundSource<float>> sources = {s};
        ffc::FreeFieldCalculator<float> calc(10, sources, c);
        ffc::FreeFieldCalculator<float>::Result r(c.height_px * c.width_px);
        calc.calculate(r);

        EXPECT_NEAR(spl,       ffc::convertMagtoSpl(r[10]), 0.1f);
        EXPECT_NEAR(spl-6.02,  ffc::convertMagtoSpl(r[20]), 0.1f);
        EXPECT_NEAR(spl-12.04, ffc::convertMagtoSpl(r[40]), 0.1f);

        auto spl_result = ffc::convert_to_spl(r);
        EXPECT_TRUE(spl_result.size() == r.size());
        for (size_t i=0; i<spl_result.size(); ++i) EXPECT_EQ(spl_result[i], ffc::convertMagtoSpl(r[i]));
    }

    TEST(LineArrayTest, setup_linearray) {
        ffc::Polar<float> p;
        LineArray<float> la {3};
        la.elementHeight = 0.1f;
        la.elementDepth = 0.25f;
        la.polar = &p;
        la.setValidSplayAngles({0.f, 5.f});
        la.recalc({0.f, 0.f});

        EXPECT_NEAR(la.elements[0].topRiggingXY.first, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[0].topRiggingXY.second, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[0].cumulativeAngle, 0.f, 1e-3f);

        EXPECT_NEAR(la.elements[1].topRiggingXY.first, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[1].topRiggingXY.second, -0.1f, 1e-3f);
        EXPECT_NEAR(la.elements[1].cumulativeAngle, 0.f, 1e-3f);

        EXPECT_NEAR(la.elements[2].topRiggingXY.first, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[2].topRiggingXY.second, -0.2f, 1e-3f);
        EXPECT_NEAR(la.elements[2].cumulativeAngle, 0.f, 1e-3f);

        la.elements[0].splayAngle = 5.f;
        la.recalc({0.f, 0.f});

        EXPECT_NEAR(la.elements[0].topRiggingXY.first, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[0].topRiggingXY.second, 0.f, 1e-3f);
        EXPECT_NEAR(la.elements[0].cumulativeAngle, 5.f, 1e-3f);

        EXPECT_NEAR(la.elements[1].bottomRiggingXY.first, -0.0087f, 1e-3f);
        EXPECT_NEAR(la.elements[1].bottomRiggingXY.second, -0.09961, 1e-3f);
        EXPECT_NEAR(la.elements[1].cumulativeAngle, 5.f, 1e-3f);
         
        // EXPECT_NEAR(geometry[2].topRiggingXY.first, 0.f, 1e-3f);
        // EXPECT_NEAR(geometry[2].topRiggingXY.second, 0.2f, 1e-3f);
    }
}