using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        public static double sinh(double x)
        {
            return (Math.Exp(x) - Math.Exp(-x)) / 2;
        }

        public static double cosh(double x)
        {
            return Math.Exp(x) / 2 + Math.Exp(-x) / 2;
        }

        public static double tanh(double x)
        {
            return sinh(x) / cosh(x);
        }

        public static double coth(double x)
        {
            return cosh(x) / sinh(x);
        }

        public static double sech(double x)
        {
            return 1.0 / cosh(x);
        }

        public static double csch(double x)
        {
            return 1.0 / sinh(x);
        }

        public static double acosh(double x)
        {
            if (x < 1)
            {
                throw new Exception("acosh requires x >= 1, but got x = " + x.ToString());
            }
            else if ((x - 1) >= XMath.epsilon)
            {
                if (x > 1 / XMath.epsilon)
                {
                    // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/06/01/06/01/0001/
                    // approximation by laurent series in 1/x at 0+ order from -1 to 0
                    return (Math.Log(x * 2));
                }
                else if (x < 1.5)
                {
                    // This is just a rearrangement of the standard form below
                    // devised to minimse loss of precision when x ~ 1:
                    double y = x - 1;
                    return log1p(y + Math.Sqrt(y * y + 2 * y));
                }
                else
                {
                    // http://functions.wolfram.com/ElementaryFunctions/ArcCosh/02/
                    return (Math.Log(x + Math.Sqrt(x * x - 1)));
                }
            }
            else
            {
                // see http://functions.wolfram.com/ElementaryFunctions/ArcCosh/06/01/04/01/0001/
                double y = x - 1;

                // approximation by taylor series in y at 0 up to order 2
                double result = Math.Sqrt(2 * y) * (1 + y / 12 + 3 * y * y / 160);
                return result;
            }
        }

        public static double asinh(double x)
        {
            if (x >= forth_root_epsilon)
            {
                if (x > 1 / root_epsilon)
                {
                    // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/06/01/06/01/0001/
                    // approximation by laurent series in 1/x at 0+ order from -1 to 1
                    return Math.Log(x * 2) + 1 / (4 * x * x);
                }
                else if (x < 0.5)
                {
                    // As below, but rearranged to preserve digits:
                    return log1p(x + sqrt1pm1(x * x));
                }
                else
                {
                    // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/02/
                    return (Math.Log(x + Math.Sqrt(x * x + 1)));
                }
            }
            else if (x <= -forth_root_epsilon)
            {
                return (-asinh(-x));
            }
            else
            {
                // http://functions.wolfram.com/ElementaryFunctions/ArcSinh/06/01/03/01/0001/
                // approximation by taylor series in x at 0 up to order 2
                double result = x;

                if (Math.Abs(x) >= root_epsilon)
                {
                    double x3 = x * x * x;
                    // approximation by taylor series in x at 0 up to order 4
                    result -= x3 / 6.0;
                }
                return (result);
            }
        }

        public static double atanh(double x)
        {
            if (x < -1) throw new Exception("atanh requires x >= -1, but got x = " + x.ToString());
            else if (x < -1 + XMath.epsilon) return double.NegativeInfinity;
            else if (x > 1 - XMath.epsilon) return double.PositiveInfinity;
            else if (x > 1) throw new Exception("atanh requires x <= 1, but got x = " + x.ToString());
            else if (Math.Abs(x) >= forth_root_epsilon)
            {
                // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/02/
                if (Math.Abs(x) < 0.5)
                    return (log1p(x) - log1p(-x)) / 2.0;
                return (Math.Log((1.0 + x) / (1.0 - x)) / 2.0);
            }
            else
            {
                // http://functions.wolfram.com/ElementaryFunctions/ArcTanh/06/01/03/01/
                // approximation by taylor series in x at 0 up to order 2
                double result = x;

                if (Math.Abs(x) >= root_epsilon)
                {
                    double x3 = x * x * x;

                    // approximation by taylor series in x at 0 up to order 4
                    result += x3 / 3.0;
                }

                return (result);
            }
        }

        public static double acoth(double x)
        {
            if (x >= -1 && x <= 1) throw new ArgumentException("Acoth is not defined for values of x between -1 and 1");
            return atanh(1 / x);
        }

        public static double asech(double x)
        {
            if (x <= 0 || x > 1) throw new ArgumentException("Asech is only defined for 0 < x <= 1");
            return acosh(1 / x);
        }

        public static double acsch(double x)
        {
            if (x == 0) throw new ArgumentException("Acsch is not defined for x = 0");
            return asinh(1 / x);
        }


        //trigonometric functions sec, cosec, acot, asec, acosec
    }
}
