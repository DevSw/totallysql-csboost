using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        public static double inverse_negative_binomial_cornish_fisher(double n, double sf, double sfc, double p, double q)
        {
            // mean:
            double m = n * (sfc) / sf;
            double t = Math.Sqrt(n * (sfc));
            // standard deviation:
            double sigma = t / sf;
            // skewness
            double sk = (1 + sfc) / t;
            // kurtosis:
            double k = (6 - sf * (5 + sfc)) / (n * (sfc));
            // Get the inverse of a std normal distribution:
            double x = erfc_inv(p > q ? 2 * q : 2 * p) * XMath.root_two;
            // Set the sign:
            if (p < 0.5)
                x = -x;
            double x2 = x * x;
            // w is correction term due to skewness
            double w = x + sk * (x2 - 1) / 6;
            //
            // Add on correction due to kurtosis.
            //
            if (n >= 10)
                w += k * x * (x2 - 3) / 24 + sk * sk * x * (2 * x2 - 5) / -36;

            w = m + sigma * w;
            if (w < XMath.min_value) return XMath.min_value;
            return w;
        }

        public static double inverse_binomial_cornish_fisher(double n, double sf, double p, double q)
        {
            // mean:
            double m = n * sf;
            // standard deviation:
            double sigma = Math.Sqrt(n * sf * (1 - sf));
            // skewness
            double sk = (1 - 2 * sf) / sigma;
            // kurtosis:
            // double  k = (1 - 6 * sf * (1 - sf) ) / (n * sf * (1 - sf));
            // Get the inverse of a std normal distribution:
            double x = XMath.erfc_inv(p > q ? 2 * q : 2 * p) * XMath.root_two;
            // Set the sign:
            if (p < 0.5)
                x = -x;
            double x2 = x * x;
            // w is correction term due to skewness
            double w = x + sk * (x2 - 1) / 6;
            /*
            // Add on correction due to kurtosis.
            // Disabled for now, seems to make things worse?
            //
            if(n >= 10)
               w += k * x * (x2 - 3) / 24 + sk * sk * x * (2 * x2 - 5) / -36;
               */
            w = m + sigma * w;
            if (w < XMath.min_value) return Math.Sqrt(XMath.min_value);
            if (w > n) return n;
            return w;
        }




    }
}
