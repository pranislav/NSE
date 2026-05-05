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
    enum Component
    {
      velocity_x = 0,
      velocity_y = 1,
      pressure   = 2
    };

    template <int dim>
    class Solution : public dealii::Function<dim>
    {
    public:
      Solution()
        : dealii::Function<dim>(dim + 1)
      {
        static_assert(dim == 2, "The MMS solution is implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < dim + 1,
                    dealii::ExcIndexRange(component, 0, dim + 1));

        if (component == velocity_x)
          return mms::velocity_x(p[0], p[1]);
        if (component == velocity_y)
          return mms::velocity_y(p[0], p[1]);
        if (component == pressure)
          return mms::pressure(p[0], p[1]);

        DEAL_II_NOT_IMPLEMENTED();
        return mms::pressure(p[0], p[1]);
      }
    };

    template <int dim>
    class RightHandSide : public dealii::Function<dim>
    {
    public:
      RightHandSide(const double viscosity)
        : dealii::Function<dim>(dim + 1)
        , viscosity(viscosity)
      {
        static_assert(dim == 2, "The MMS right hand side is implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < dim + 1,
                    dealii::ExcIndexRange(component, 0, dim + 1));

        if (component == velocity_x)
          return mms::momentum_rhs_x(p[0], p[1], viscosity);
        if (component == velocity_y)
          return mms::momentum_rhs_y(p[0], p[1], viscosity);
        if (component == pressure)
          return mms::continuity_rhs(p[0], p[1]);

        DEAL_II_NOT_IMPLEMENTED();
        return mms::continuity_rhs(p[0], p[1]);
      }

    private:
      const double viscosity;
    };

    template <int dim>
    class TemperatureSolution : public dealii::Function<dim>
    {
    public:
      TemperatureSolution()
        : dealii::Function<dim>(1)
      {
        static_assert(dim == 2,
                      "The MMS temperature solution is implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature(p[0], p[1]);
      }
    };

    template <int dim>
    class TemperatureRightHandSide : public dealii::Function<dim>
    {
    public:
      TemperatureRightHandSide(const double thermal_diffusivity)
        : dealii::Function<dim>(1)
        , thermal_diffusivity(thermal_diffusivity)
      {
        static_assert(dim == 2,
                      "The MMS temperature right hand side is implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature_rhs(p[0], p[1], thermal_diffusivity);
      }

    private:
      const double thermal_diffusivity;
    };

    template <int dim>
    class VelocityBoundaryValues : public dealii::Function<dim>
    {
    public:
      VelocityBoundaryValues()
        : dealii::Function<dim>(dim + 1)
      {
        static_assert(dim == 2,
                      "The MMS velocity boundary values are implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component) const override
      {
        AssertThrow(component < dim + 1,
                    dealii::ExcIndexRange(component, 0, dim + 1));

        if (component == velocity_x)
          return mms::velocity_x(p[0], p[1]);
        if (component == velocity_y)
          return mms::velocity_y(p[0], p[1]);
        if (component == pressure)
          return 0.0;

        DEAL_II_NOT_IMPLEMENTED();
        return 0.0;
      }
    };

    template <int dim>
    class TemperatureBoundaryValues : public dealii::Function<dim>
    {
    public:
      TemperatureBoundaryValues()
        : dealii::Function<dim>(1)
      {
        static_assert(dim == 2,
                      "The MMS temperature boundary values are implemented only in 2D.");
      }

      virtual double value(const dealii::Point<dim> &p,
                           const unsigned int        component = 0) const override
      {
        AssertThrow(component < 1, dealii::ExcIndexRange(component, 0, 1));
        return mms::temperature(p[0], p[1]);
      }
    };
  } // namespace MMS
} // namespace Cht

#endif
