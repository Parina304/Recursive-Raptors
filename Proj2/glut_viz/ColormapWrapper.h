#include "cppcolormap.h"
#include <xtensor/xview.hpp>
#include <stdexcept>

/**
 * Get the color from a colormap at a specific position.
 * 
 * @param colormap The colormap function (e.g., cppcolormap::jet).
 * @param position A value between 0 and 1 representing the position in the colormap.
 * @param steps The number of steps to interpolate the colormap (default: 256).
 * @return A 1D xtensor containing the RGB values of the color.
 */

 

template <class T, class R = xt::xtensor<double, 2>>
inline R interp(const T& arg, size_t N)
{
    CPPCOLORMAP_ASSERT(arg.dimension() == 2);
    using size_type = typename T::shape_type::value_type;
    size_type n = static_cast<size_type>(N);
    size_type m = static_cast<size_type>(arg.shape(1));

    if (n == arg.shape(0)) {
        return arg;
    }

    R ret({n, m});

    xt::xtensor<double, 1> x = xt::linspace(0.0, 1.0, arg.shape(0));
    xt::xtensor<double, 1> xi = xt::linspace(0.0, 1.0, n);

    for (size_t j = 0; j < arg.shape(1); j++) {
        auto c = xt::view(arg, xt::all(), j);
        auto ci = xt::view(ret, xt::all(), j);
        ci = xt::interp(xi, x, c);
    }

    return ret;
}

inline xt::xtensor<double, 1> get_color_from_colormap(
    const std::function<xt::xtensor<double, 2>(size_t)>& colormap,
    double position,
    size_t steps = 256)
{
    if (position < 0.0 || position > 1.0) {
        throw std::out_of_range("Position must be between 0 and 1.");
    }

    // Generate the colormap with the specified number of steps
    auto cmap = colormap(steps);

    // Calculate the index corresponding to the position
    size_t index = static_cast<size_t>(position * (cmap.shape(0) - 1));

    // Extract and return the color at the calculated index
    return xt::view(cmap, index, xt::all());
}

