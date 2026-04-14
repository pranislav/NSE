#ifndef STEP57_BOUNDARY_CONDITIONS_H
#define STEP57_BOUNDARY_CONDITIONS_H

#include "case_config.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

namespace Cht
{
  template <int dim>
  class VelocityBoundaryValues : public dealii::Function<dim>
  {
  public:
    VelocityBoundaryValues(const VelocityBoundaryCondition &condition,
                           const double                    coordinate_min,
                           const double                    coordinate_max)
      : dealii::Function<dim>(dim + 1)
      , condition(condition)
      , coordinate_min(coordinate_min)
      , coordinate_max(coordinate_max)
    {}

    virtual double value(const dealii::Point<dim> &p,
                         const unsigned int        component) const override
    {
      if (component != condition.component)
        return 0.0;

      if (condition.type == VelocityBoundaryCondition::Type::constant)
        return condition.value;

      const double width = coordinate_max - coordinate_min;
      AssertThrow(width > 0.0,
                  dealii::ExcMessage(
                    "Parabolic velocity boundary has zero width."));

      const double s = p[condition.coordinate];
      if (s < coordinate_min || s > coordinate_max)
        return 0.0;

      const double normalized = (s - coordinate_min) * (coordinate_max - s);
      return 4.0 * condition.value * normalized / (width * width);
    }

  private:
    const VelocityBoundaryCondition condition;
    const double                    coordinate_min;
    const double                    coordinate_max;
  };

  template <int dim>
  class TemperatureBoundaryValues : public dealii::Function<dim>
  {
  public:
    explicit TemperatureBoundaryValues(const double temperature)
      : dealii::Function<dim>(1)
      , temperature(temperature)
    {}

    virtual double value(const dealii::Point<dim> &,
                         const unsigned int        component) const override
    {
      AssertIndexRange(component, 1);
      return temperature;
    }

  private:
    const double temperature;
  };
}

#endif
