#ifndef STEP57_MMS_H
#define STEP57_MMS_H

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

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

      virtual dealii::Tensor<1, dim> gradient(const dealii::Point<dim> &p,
                              const unsigned int        component = 0) const override
        {
          AssertThrow(component < dim + 1,
                      dealii::ExcIndexRange(component, 0, dim + 1));

          dealii::Tensor<1, 2> return_value;
          
          if (component == velocity_x)
            return_value[0] = mms::grad_x_velocity_x(p[0], p[1]);
            return_value[1] = mms::grad_y_velocity_x(p[0], p[1]);
          if (component == velocity_y)
            return_value[0] = mms::grad_x_velocity_y(p[0], p[1]);
            return_value[1] = mms::grad_y_velocity_x(p[0], p[1]);
          if (component == pressure)
            return_value[0] = mms::grad_x_pressure(p[0], p[1]);
            return_value[1] = mms::grad_y_pressure(p[0], p[1]);

          return return_value;
        }
    };

    template <int dim>
    class RightHandSide
    {
    public:
      RightHandSide(const double viscosity)
        : viscosity(viscosity)
      {
        static_assert(dim == 2, "The MMS right hand side is implemented only in 2D.");
      }

      dealii::Tensor<1, dim> value(const dealii::Point<dim> &p) const
      {
        dealii::Tensor<1, dim> rhs;
        rhs[velocity_x] = mms::momentum_rhs_x(p[0], p[1], viscosity);
        rhs[velocity_y] = mms::momentum_rhs_y(p[0], p[1], viscosity);
        return rhs;
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
          return 0.0; // assertion was not passing thorugh dealii vector_value

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
