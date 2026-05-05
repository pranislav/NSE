#ifndef STEP57_MMS_H
#define STEP57_MMS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include "mms_generated.h"

namespace Cht
{
  namespace MMS
  {
    template <int dim>
    class Solution;

    template <>
    class Solution<2> : public dealii::Function<2>
    {
    public:
      Solution()
        : dealii::Function<2>(2 + 1)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 2 + 1, dealii::ExcIndexRange(component, 0, 2 + 1));

        if (component == 0)
          return mms::velocity_x(p[0], p[1]);
        if (component == 1)
          return mms::velocity_y(p[0], p[1]);
        return mms::pressure(p[0], p[1]);
      }
    };

    template <int dim>
    class RightHandSide;

    template <>
    class RightHandSide<2> : public dealii::Function<2>
    {
    public:
      RightHandSide(const double viscosity)
        : dealii::Function<2>(2 + 1)
        , viscosity(viscosity)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 2 + 1, dealii::ExcIndexRange(component, 0, 2 + 1));

        if (component == 0)
          return mms::momentum_rhs_x(p[0], p[1], viscosity);
        if (component == 1)
          return mms::momentum_rhs_y(p[0], p[1], viscosity);
        return mms::continuity_rhs(p[0], p[1]);
      }

    private:
      const double viscosity;
    };

    template <int dim>
    class TemperatureSolution;

    template <>
    class TemperatureSolution<2> : public dealii::Function<2>
    {
    public:
      TemperatureSolution()
        : dealii::Function<2>(1)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature(p[0], p[1]);
      }
    };

    template <int dim>
    class TemperatureRightHandSide;

    template <>
    class TemperatureRightHandSide<2> : public dealii::Function<2>
    {
    public:
      TemperatureRightHandSide(const double thermal_diffusivity)
        : dealii::Function<2>(1)
        , thermal_diffusivity(thermal_diffusivity)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature_rhs(p[0], p[1], thermal_diffusivity);
      }

    private:
      const double thermal_diffusivity;
    };

    template <int dim>
    class VelocityBoundaryValues;

    template <>
    class VelocityBoundaryValues<2> : public dealii::Function<2>
    {
    public:
      VelocityBoundaryValues()
        : dealii::Function<2>(2 + 1)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component) const override
      {
        AssertThrow(component < 2 + 1, dealii::ExcIndexRange(component, 0, 2 + 1));

        if (component == 0)
          return mms::velocity_x(p[0], p[1]);
        if (component == 1)
          return mms::velocity_y(p[0], p[1]);
        return 0.0;
      }
    };

    template <int dim>
    class TemperatureBoundaryValues;

    template <>
    class TemperatureBoundaryValues<2> : public dealii::Function<2>
    {
    public:
      TemperatureBoundaryValues()
        : dealii::Function<2>(1)
      {}

      virtual double value(const dealii::Point<2> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature(p[0], p[1]);
      }
    };
  } // namespace MMS
} // namespace Cht

#endif
