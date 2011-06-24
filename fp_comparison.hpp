//  boost fp_comparison.hpp header file

// (C) Copyright Alberto Ganesh Barbati 2006.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)

// See http://www.boost.org/libs/math for the library home page.
// See http://www.boost.org for updates, documentation, and revision history.

#ifndef BOOST_MATH_FP_COMPARISON_HPP
#define BOOST_MATH_FP_COMPARISON_HPP

#include <cmath>
#include <limits>
#include <cassert>
#include <functional>

#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace boost
{
    namespace math
    {
        template <typename T = double>
        struct euclidean_distance : std::binary_function<T, T, T>
        {
            T operator()(T x, T y) const
            {
                using std::abs;
                return abs(x - y);
            }
        };

        namespace detail
        {
            // copied from <boost/test/floating_point_comparison.hpp>
            template <typename T>
            inline T safe_fp_division(T x, T y)
            {
                // both x and y are unsigned here
                return y < 1 && x > y * (std::numeric_limits<T>::max)()             ? (std::numeric_limits<T>::max)()
                     : y > 1 && x < y * (std::numeric_limits<T>::min)() || x == 0   ? 0
                     :                                                                x / y;
            }
        }

        template <typename T = double>
        struct relative_error : std::binary_function<T, T, T>
        {
            T operator()(T x, T y) const
            {
                using std::abs;
                return detail::safe_fp_division(abs(x - y), abs(x));
            }
        };

        template <typename T = double>
        struct relative_error_fast : std::binary_function<T, T, T>
        {
            T operator()(T x, T y) const
            {
                using std::abs;
                return abs((x - y) / x);
            }
        };

        template <class Distance>
        struct symmetric_distance_adaptor : private Distance
        {
            // behaves like std::binary_function
            typedef BOOST_DEDUCED_TYPENAME Distance::first_argument_type    first_argument_type;
            typedef BOOST_DEDUCED_TYPENAME Distance::second_argument_type   second_argument_type;
            typedef BOOST_DEDUCED_TYPENAME Distance::result_type            result_type;

            BOOST_STATIC_ASSERT((boost::is_same<first_argument_type, second_argument_type>::value));

            symmetric_distance_adaptor()
            {}

            symmetric_distance_adaptor(Distance d)
                : Distance(d)
            {}

            result_type operator()(first_argument_type x, second_argument_type y) const
            {
                return (std::max)(Distance::operator()(x, y), Distance::operator()(y, x));
            }
        };

        template <typename T = double>
        struct relative_error_symmetric : symmetric_distance_adaptor<relative_error<T> >
        {};

        template <typename T = double>
        struct relative_error_symmetric_fast : symmetric_distance_adaptor<relative_error_fast<T> >
        {};

        template <typename T = double>
        struct s1_distance : std::binary_function<T, T, T>
        {
            T period_;

            s1_distance(T period)
                : period_(period)
            {
                assert(period_ > 0);
            }

            T operator()(T x, T y) const
            {
                using std::abs;
                assert(abs(x - y) <= period_);
                return (x < y)
                    ? (std::min)(y - x, x - y + period_)
                    : (std::min)(x - y, y - x + period_);
            }
        };

        template <class Distance = euclidean_distance<double> >
        struct are_close : private Distance
        {
            // behaves like std::binary_function
            typedef BOOST_DEDUCED_TYPENAME Distance::first_argument_type    first_argument_type;
            typedef BOOST_DEDUCED_TYPENAME Distance::second_argument_type   second_argument_type;
            typedef bool                                                    result_type;

            typedef BOOST_DEDUCED_TYPENAME Distance::result_type            distance_type;

            distance_type epsilon_;

            are_close(distance_type epsilon)
                : Distance(), epsilon_(epsilon)
            {
                assert(epsilon_ > 0);
            }

            template <typename T1>
            are_close(distance_type epsilon, const T1& t1)
                : Distance(t1), epsilon_(epsilon)
            {}

            template <typename T1, typename T2>
            are_close(distance_type epsilon, const T1& t1, const T2& t2)
                : Distance(t1, t2), epsilon_(epsilon)
            {}

            template <typename T1, typename T2, typename T3>
            are_close(distance_type epsilon, const T1& t1, const T2& t2, const T3& t3)
                : Distance(t1, t2, t3), epsilon_(epsilon)
            {}

            bool operator()(first_argument_type x, second_argument_type y) const
            {
                return Distance::operator()(x, y) < epsilon_;
            }
        };

        // cheap substitute to std::bind1st(are_close<>(), origin)

        template <class Distance = euclidean_distance<double> >
        struct is_close_to : private Distance
        {
            // behaves like std::unary_function
            typedef BOOST_DEDUCED_TYPENAME Distance::second_argument_type   argument_type;
            typedef bool                                                    result_type;

            typedef BOOST_DEDUCED_TYPENAME Distance::first_argument_type    origin_type;
            typedef BOOST_DEDUCED_TYPENAME Distance::result_type            distance_type;

            origin_type     origin_;
            distance_type   epsilon_;

            is_close_to(origin_type origin, distance_type epsilon)
                : Distance(), origin_(origin), epsilon_(epsilon)
            {}

            template <typename T1>
            is_close_to(origin_type origin, distance_type epsilon, const T1& t1)
                : Distance(t1), origin_(origin), epsilon_(epsilon)
            {}

            template <typename T1, typename T2>
            is_close_to(origin_type origin, distance_type epsilon, const T1& t1, const T2& t2)
                : Distance(t1, t2), origin_(origin), epsilon_(epsilon)
            {}

            template <typename T1, typename T2, typename T3>
            is_close_to(origin_type origin, distance_type epsilon, const T1& t1, const T2& t2, const T3& t3)
                : Distance(t1, t2, t3), origin_(origin), epsilon_(epsilon)
            {}

            bool operator()(argument_type x) const
            {
                return Distance::operator()(origin_, x) < epsilon_;
            }
        };

    } // namespace math

} // namespace boost

#endif // BOOST_MATH_FP_COMPARISON_HPP