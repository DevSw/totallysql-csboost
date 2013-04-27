using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    partial class XMath
    {
        public static double inverse_poisson_cornish_fisher(double lambda, double p, double q)
        {
            double m = lambda;
            // standard deviation:
            double sigma = Math.Sqrt(lambda);
            // skewness
            double sk = 1 / sigma;
            // kurtosis:
            // double k = 1/lambda;
            // Get the inverse of a std normal distribution:
            double x = erfc_inv(p > q ? 2 * q : 2 * p) * XMath.root_two;
            // Set the sign:
            if (p < 0.5) x = -x;
            double x2 = x * x;
            // w is correction term due to skewness
            double w = x + sk * (x2 - 1) / 6;
            w = m + sigma * w;
            return w > XMath.min_value ? w : XMath.min_value;
        }

    }
}
