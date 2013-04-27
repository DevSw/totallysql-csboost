using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        const int pn_Size = 30;

        public static double beta(double a, double b)
        {
            if (a <= 0)
                throw new Exception(string.Format("The arguments to the beta function must be greater than zero (got a={0:G}).", a));
            if (b <= 0)
                throw new Exception(string.Format("The arguments to the beta function must be greater than zero (got b={0:G}).", b));

            double result;

            double c = a + b;

            // special cases:
            if ((c == a) && (b < XMath.epsilon)) return gamma(b);
            else if ((c == b) && (a < XMath.epsilon)) return gamma(a);
            else if (b == 1) return 1 / a;
            else if (a == 1) return 1 / b;


            if (a < b)
            {
                double d = a;
                a = b;
                b = d;
            }

            double agh = a + lanczos_g - 0.5;
            double bgh = b + lanczos_g - 0.5;
            double cgh = c + lanczos_g - 0.5;
            result = lanczos_sum_expG_scaled(a) * lanczos_sum_expG_scaled(b) / lanczos_sum_expG_scaled(c);
            double ambh = a - 0.5 - b;

            if ((Math.Abs(b * ambh) < (cgh * 100)) && (a > 100))
            {
                // Special case where the base of the power term is close to 1
                // compute (1+x)^y instead:
                result *= Math.Exp(ambh * log1p(-b / cgh));
            }
            else
            {
                result *= Math.Pow(agh / cgh, a - 0.5 - b);
            }
            if (cgh > 1e10)
                // this avoids possible overflow, but appears to be marginally less accurate:
                result *= Math.Pow((agh / cgh) * (bgh / cgh), b);
            else
                result *= Math.Pow((agh * bgh) / (cgh * cgh), b);
            result *= Math.Sqrt(Math.E / bgh);

            return result;
        }

        public static double ibetac(double a, double b, double x)
        {
            double p_derivative = 0.0;
            return ibeta_imp(a, b, x, true, true, ref p_derivative);
        }

        public static double betac(double a, double b, double x)
        {
            double p_derivative = 0.0;
            return ibeta_imp(a, b, x, true, false, ref p_derivative);
        }

        public static double ibeta(double a, double b, double x)
        {
            double p_derivative = 0.0;
            return ibeta_imp(a, b, x, false, true, ref p_derivative);
        }

        public static double beta(double a, double b, double x)
        {
            double p_derivative = 0.0;
            return ibeta_imp(a, b, x, false, false, ref p_derivative);
        }

        public static double ibeta_derivative(double a, double b, double x)
        {
            if (a <= 0)
                throw new Exception(string.Format("The argument a to the incomplete beta function must be greater than zero (got a={0:G}).", a));
            if (b <= 0)
                throw new Exception(string.Format("The argument b to the incomplete beta function must be greater than zero (got b={0:G}).", b));
            if ((x < 0) || (x > 1))
                throw new Exception(string.Format("Parameter x outside the range [0,1] in the incomplete beta function (got x={0:G}).", x));
            //
            // Now the corner cases:
            //
            if (x == 0)
            {
                if (a > 1) return 0.0;
                if (a == 1) return 1.0 / beta(a, b);
                throw new OverflowException();
            }
            else if (x == 1)
            {
                if (b > 1) return 0.0;
                if (b == 1) return 1.0 / beta(a, b);
                throw new OverflowException();
            }
            //
            // Now the regular cases:
            //
            double f1 = ibeta_power_terms(a, b, x, 1 - x, true);
            double y = (1 - x) * x;

            if (f1 == 0) return 0;

            if (double.MaxValue * y < f1) throw new OverflowException();

            f1 /= y;

            return f1;
        }

        public static double ibeta_inv(double a, double b, double p)
        {
            if (a <= 0)
                throw new Exception(string.Format("The argument a to the incomplete beta function inverse must be greater than zero (got a={0:G}).", a));
            if (b <= 0)
                throw new Exception(string.Format("The argument b to the incomplete beta function inverse must be greater than zero (got b={0:G}).", b));
            if ((p < 0) || (p > 1))
                throw new Exception(string.Format("Argument p outside the range [0,1] in the incomplete beta function inverse (got p={0:G}).", p));

            double rx;
            double py = 0;
            rx = ibeta_inv_imp(a, b, p, 1.0 - p, ref py);
            return rx;
        }

        public static double ibeta_invb(double a, double x, double p)
        {
            if (p == 0) return XMath.min_value;
            if (p == 1) return double.MaxValue;

            return ibeta_inv_ab_imp(a, x, p, 1 - p, true);
        }

        public static double ibeta_inva(double b, double x, double p)
        {
            if (p == 0) return double.MaxValue;
            if (p == 1) return XMath.min_value;

            return ibeta_inv_ab_imp(b, x, p, 1 - p, false);
        }

        public static double ibetac_inv(double a, double b, double q)
        {
            if (a <= 0)
                throw new Exception(string.Format("The argument a to the incomplete beta function inverse must be greater than zero (got a={0:G}).", a));
            if (b <= 0)
                throw new Exception(string.Format("The argument b to the incomplete beta function inverse must be greater than zero (got b={0:G}).", b));
            if ((q < 0) || (q > 1))
                throw new Exception(string.Format("Argument q outside the range [0,1] in the incomplete beta function inverse (got p={0:G}).", q));

            double rx;
            double py = 0;

            rx = ibeta_inv_imp(a, b, 1.0 - q, q, ref py);
            return rx;
        }

        public static double ibetac_invb(double a, double x, double q)
        {
            if (q == 1) return XMath.min_value;
            if (q == 0) return double.MaxValue;

            return ibeta_inv_ab_imp(a, x, 1 - q, q, true);
        }

        public static double ibetac_inva(double b, double x, double q)
        {
            if (q == 1) return double.MaxValue;
            if (q == 0) return XMath.min_value;

            return ibeta_inv_ab_imp(b, x, 1 - q, q, false);
        }

        #region helper functions

        private class ibeta_series_t : series<double>
        {
            double result, x, apn, poch;
            int n;

            public ibeta_series_t(double a_, double b_, double x_, double mult)
            {
                result = mult;
                x = x_;
                apn = a_;
                poch = 1 - b_;
                n = 1;
            }
            public override double next()
            {
                double r = result / apn;
                apn += 1;
                result *= poch * x / n;
                ++n;
                poch += 1;
                return r;
            }
        }

        private static double beta_small_b_large_a_series(double a, double b, double x, double y, double s0, double mult, bool normalised)
        {
            //
            // This is DiDonato and Morris's BGRAT routine, see Eq's 9 through 9.6.
            //
            // Some values we'll need later, these are Eq 9.1:
            //
            double bm1 = b - 1;
            double t = a + bm1 / 2;
            double lx, u;
            if (y < 0.35) lx = log1p(-y);
            else lx = Math.Log(x);
            u = -t * lx;
            // and from from 9.2:
            double prefix;
            double h = regularised_gamma_prefix(b, u);
            if (h <= XMath.min_value) return s0;
            if (normalised)
            {
                prefix = h / gamma_delta_ratio(a, b);
                prefix /= Math.Pow(t, b);
            }
            else
            {
                prefix = full_igamma_prefix(b, u) / Math.Pow(t, b);
            }
            prefix *= mult;
            //
            // now we need the quantity Pn, unfortunatately this is computed
            // recursively, and requires a full history of all the previous values
            // so no choice but to declare a big table and hope it's big enough...
            //
            double[] p = new double[pn_Size];  // see 9.3.
            p[0] = 1;
            //
            // Now an initial value for J, see 9.6:
            //
            double j = gamma_q(b, u) / h;
            //
            // Now we can start to pull things together and evaluate the sum in Eq 9:
            //
            double sum = s0 + prefix * j;  // Value at N = 0
            // some variables we'll need:
            uint tnp1 = 1; // 2*N+1
            double lx2 = lx / 2;
            lx2 *= lx2;
            double lxp = 1;
            double t4 = 4 * t * t;
            double b2n = b;

            for (uint n = 1; n < p.Length; ++n)
            {
                /*
                // debugging code, enable this if you want to determine whether
                // the table of Pn's is large enough...
                //
                static int max_count = 2;
                if(n > max_count)
                {
                   max_count = n;
                   std::cerr << "Max iterations in BGRAT was " << n << std::endl;
                }
                */
                //
                // begin by evaluating the next Pn from Eq 9.4:
                //
                tnp1 += 2;
                p[n] = 0;
                double mbn = b - n;
                uint tmp1 = 3;
                double[] f = factorials();
                for (uint m = 1; m < n; ++m)
                {
                    mbn = m * b - n;
                    p[n] += mbn * p[n - m] / f[tmp1];
                    tmp1 += 2;
                }
                p[n] /= n;
                p[n] += bm1 / f[tnp1];
                //
                // Now we want Jn from Jn-1 using Eq 9.6:
                //
                j = (b2n * (b2n + 1) * j + (u + b2n + 1) * lxp) / t4;
                lxp *= lx2;
                b2n += 2;
                //
                // pull it together with Eq 9:
                //
                double r = prefix * p[n] * j;
                sum += r;
                if (r > 1)
                {
                    if (Math.Abs(r) < Math.Abs(XMath.epsilon * sum))
                        break;
                }
                else
                {
                    if (Math.Abs(r / XMath.epsilon) < Math.Abs(sum))
                        break;
                }
            }
            return sum;
        }

        private static double binomial_ccdf(double n, double k, double x, double y)
        {
            double result = Math.Pow(x, n);
            double term = result;
            for (uint i = (uint)(n - 1.0); i > k; --i)
            {
                term *= ((i + 1) * y) / ((n - i) * x);
                result += term;
            }

            return result;
        }

        private static double ibeta_power_terms(double a, double b, double x, double y, bool normalised)
        {

            if (!normalised) return Math.Pow(x, a) * Math.Pow(y, b);

            double result;

            double prefix = 1;
            double c = a + b;

            // combine power terms with Lanczos approximation:
            double agh = a + lanczos_g - 0.5;
            double bgh = b + lanczos_g - 0.5;
            double cgh = c + lanczos_g - 0.5;
            result = lanczos_sum_expG_scaled(c) / (lanczos_sum_expG_scaled(a) * lanczos_sum_expG_scaled(b));

            // l1 and l2 are the base of the exponents minus one:
            double l1 = (x * b - y * agh) / agh;
            double l2 = (y * a - x * bgh) / bgh;
            if ((min(Math.Abs(l1), Math.Abs(l2)) < 0.2))
            {
                // when the base of the exponent is very near 1 we get really
                // gross errors unless extra care is taken:
                if ((l1 * l2 > 0) || (min(a, b) < 1))
                {
                    //
                    // This first branch handles the simple cases where either: 
                    //
                    // * The two power terms both go in the same direction 
                    // (towards zero or towards infinity).  In this case if either 
                    // term overflows or underflows, then the product of the two must 
                    // do so also.  
                    // *Alternatively if one exponent is less than one, then we 
                    // can't productively use it to eliminate overflow or underflow 
                    // from the other term.  Problems with spurious overflow/underflow 
                    // can't be ruled out in this case, but it is *very* unlikely 
                    // since one of the power terms will evaluate to a number close to 1.
                    //
                    if (Math.Abs(l1) < 0.1) result *= Math.Exp(a * log1p(l1));
                    else result *= Math.Pow((x * cgh) / agh, a);
                    if (Math.Abs(l2) < 0.1) result *= Math.Exp(b * log1p(l2));
                    else result *= Math.Pow((y * cgh) / bgh, b);
                }
                else if (max(Math.Abs(l1), Math.Abs(l2)) < 0.5)
                {
                    //
                    // Both exponents are near one and both the exponents are 
                    // greater than one and further these two 
                    // power terms tend in opposite directions (one towards zero, 
                    // the other towards infinity), so we have to combine the terms 
                    // to avoid any risk of overflow or underflow.
                    //
                    // We do this by moving one power term inside the other, we have:
                    //
                    //    (1 + l1)^a * (1 + l2)^b
                    //  = ((1 + l1)*(1 + l2)^(b/a))^a
                    //  = (1 + l1 + l3 + l1*l3)^a   ;  l3 = (1 + l2)^(b/a) - 1
                    //                                    = Math.Exp((b/a) * Math.Log(1 + l2)) - 1
                    //
                    // The tricky bit is deciding which term to move inside :-)
                    // By preference we move the larger term inside, so that the
                    // size of the largest exponent is reduced.  However, that can
                    // only be done as long as l3 (see above) is also small.
                    //
                    bool small_a = a < b;
                    double ratio = b / a;
                    if ((small_a && (ratio * l2 < 0.1)) || (!small_a && (l1 / ratio > 0.1)))
                    {
                        double l3 = expm1(ratio * log1p(l2));
                        l3 = l1 + l3 + l3 * l1;
                        l3 = a * log1p(l3);
                        result *= Math.Exp(l3);
                    }
                    else
                    {
                        double l3 = expm1(log1p(l1) / ratio);
                        l3 = l2 + l3 + l3 * l2;
                        l3 = b * log1p(l3);
                        result *= Math.Exp(l3);
                    }
                }
                else if (Math.Abs(l1) < Math.Abs(l2))
                {
                    // First base near 1 only:
                    double l = a * log1p(l1)
                       + b * Math.Log((y * cgh) / bgh);
                    result *= Math.Exp(l);
                }
                else
                {
                    // Second base near 1 only:
                    double l = b * log1p(l2)
                       + a * Math.Log((x * cgh) / agh);
                    result *= Math.Exp(l);
                }
            }
            else
            {
                // general case:
                double b1 = (x * cgh) / agh;
                double b2 = (y * cgh) / bgh;
                l1 = a * Math.Log(b1);
                l2 = b * Math.Log(b2);
                if ((l1 >= log_max_value)
                   || (l1 <= log_min_value)
                   || (l2 >= log_max_value)
                   || (l2 <= log_min_value)
                   )
                {
                    // Oops, overflow, sidestep:
                    if (a < b)
                        result *= Math.Pow(Math.Pow(b2, b / a) * b1, a);
                    else
                        result *= Math.Pow(Math.Pow(b1, a / b) * b2, b);
                }
                else
                {
                    // finally the normal case:
                    result *= Math.Pow(b1, a) * Math.Pow(b2, b);
                }
            }
            // combine with the leftover terms from the Lanczos approximation:
            result *= Math.Sqrt(bgh / Math.E);
            result *= Math.Sqrt(agh / cgh);
            result *= prefix;

            return result;
        }

        private static double ibeta_series(double a, double b, double x, double s0, bool normalised, ref double p_derivative, double y)
        {

            double result;

            if (normalised)
            {
                double c = a + b;

                // incomplete beta power term, combined with the Lanczos approximation:
                double agh = a + lanczos_g - 0.5;
                double bgh = b + lanczos_g - 0.5;
                double cgh = c + lanczos_g - 0.5;
                result = lanczos_sum_expG_scaled(c) / (lanczos_sum_expG_scaled(a) * lanczos_sum_expG_scaled(b));
                if (a * b < bgh * 10)
                    result *= Math.Exp((b - 0.5) * log1p(a / bgh));
                else
                    result *= Math.Pow(cgh / bgh, b - 0.5);
                result *= Math.Pow(x * cgh / agh, a);
                result *= Math.Sqrt(agh / Math.E);

                if (p_derivative != 0) p_derivative = result * Math.Pow(y, b);
            }
            else
            {
                // Non-normalised, just compute the power:
                result = Math.Pow(x, a);
            }
            if (result < XMath.min_value) return s0; // Safeguard: series can't cope with denorms.
            ibeta_series_t s = new ibeta_series_t(a, b, x, result);
            int max_iter = max_series_iterations;
            result = sum_series(s, XMath.epsilon, ref max_iter, s0);
            return result;
        }

        private static double ibeta_a_step(double a, double b, double x, double y, int k, bool normalised, ref double p_derivative)
        {
            double prefix = ibeta_power_terms(a, b, x, y, normalised);
            if (p_derivative != 0) p_derivative = prefix;
            prefix /= a;
            if (prefix == 0) return prefix;
            double sum = 1.0;
            double term = 1.0;
            // series summation from 0 to k-1:
            for (int i = 0; i < k - 1; ++i)
            {
                term *= (a + b + i) * x / (a + i + 1);
                sum += term;
            }
            prefix *= sum;

            return prefix;
        }

        private static double rising_factorial_ratio(double a, double b, int k)
        {
            // calculate:
            // (a)(a+1)(a+2)...(a+k-1)
            // _______________________
            // (b)(b+1)(b+2)...(b+k-1)

            // This is only called with small k, for large k
            // it is grossly inefficient, do not use outside it's
            // intended purpose!!!
            if (k == 0)
                return 1;
            double result = 1.0;
            for (int i = 0; i < k; ++i)
                result *= (a + i) / (b + i);
            return result;
        }

        internal static double ibeta_imp(double a, double b, double x, bool inv, bool normalised, ref double p_derivative)
        {
            bool invert = inv;
            double fract;
            double y = 1 - x;

            p_derivative = -1; // value not set.

            if ((x < 0) || (x > 1))
                throw new Exception(string.Format("Parameter x outside the range [0,1] in the incomplete beta function (got x={0:G}).", x));

            if (normalised)
            {
                if (a < 0)
                    throw new Exception(string.Format("The argument a to the incomplete beta function must be >= zero (got a={0:G}).", a));
                if (b < 0)
                    throw new Exception(string.Format("The argument b to the incomplete beta function must be >= zero (got b={0:G}).", b));
                // extend to a few very special cases:
                if (a == 0)
                {
                    if (b == 0)
                        throw new Exception(string.Format("The arguments a and b to the incomplete beta function cannot both be zero, with x={0:G}.", x));
                    if (b > 0)
                        return inv ? 0 : 1;
                }
                else if (b == 0)
                {
                    if (a > 0)
                        return inv ? 1 : 0;
                }
            }
            else
            {
                if (a <= 0)
                    throw new Exception(string.Format("The argument a to the incomplete beta function must be greater than zero (got a={0:G}).", a));
                if (b <= 0)
                    throw new Exception(string.Format("The argument b to the incomplete beta function must be greater than zero (got b={0:G}).", b));
            }

            if (x == 0)
            {
                p_derivative = (a == 1) ? (double)1 : (a < 1) ? double.MaxValue / 2.0 : XMath.min_value * 2.0;
                return (invert ? (normalised ? 1.0 : beta(a, b)) : 0.0);
            }
            if (x == 1)
            {
                p_derivative = (b == 1) ? 1.0 : (b < 1) ? double.MaxValue / 2.0 : XMath.min_value * 2.0;
                return (!invert ? (normalised ? 1.0 : beta(a, b)) : 0.0);
            }

            if (min(a, b) <= 1)
            {
                if (x > 0.5)
                {
                    swap(ref a, ref b);
                    swap(ref x, ref y);
                    invert = !invert;
                }
                if (max(a, b) <= 1)
                {
                    // Both a,b < 1:
                    if ((a >= min(0.2, b)) || (Math.Pow(x, a) <= 0.9))
                    {
                        if (!invert)
                        {
                            fract = ibeta_series(a, b, x, 0.0, normalised, ref p_derivative, y);
                        }
                        else
                        {
                            fract = -(normalised ? 1 : beta(a, b));
                            invert = false;
                            fract = -ibeta_series(a, b, x, fract, normalised, ref p_derivative, y);
                        }
                    }
                    else
                    {
                        swap(ref a, ref b);
                        swap(ref x, ref y);
                        invert = !invert;
                        if (y >= 0.3)
                        {
                            if (!invert)
                            {
                                fract = ibeta_series(a, b, x, 0.0, normalised, ref p_derivative, y);
                            }
                            else
                            {
                                fract = -(normalised ? 1 : beta(a, b));
                                invert = false;
                                fract = -ibeta_series(a, b, x, fract, normalised, ref p_derivative, y);
                            }
                        }
                        else
                        {
                            // Sidestep on a, and then use the series representation:
                            double prefix;
                            if (!normalised)
                            {
                                prefix = rising_factorial_ratio(a + b, a, 20);
                            }
                            else
                            {
                                prefix = 1;
                            }
                            fract = ibeta_a_step(a, b, x, y, 20, normalised, ref p_derivative);
                            if (!invert)
                            {
                                fract = beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, normalised);
                            }
                            else
                            {
                                fract -= (normalised ? 1 : beta(a, b));
                                invert = false;
                                fract = -beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, normalised);
                            }
                        }
                    }
                }
                else
                {
                    // One of a, b < 1 only:
                    if ((b <= 1) || ((x < 0.1) && (Math.Pow(b * x, a) <= 0.7)))
                    {
                        if (!invert)
                        {
                            fract = ibeta_series(a, b, x, 0.0, normalised, ref p_derivative, y);
                        }
                        else
                        {
                            fract = -(normalised ? 1 : beta(a, b));
                            invert = false;
                            fract = -ibeta_series(a, b, x, fract, normalised, ref p_derivative, y);
                        }
                    }
                    else
                    {
                        swap(ref a, ref b);
                        swap(ref x, ref y);
                        invert = !invert;

                        if (y >= 0.3)
                        {
                            if (!invert)
                            {
                                fract = ibeta_series(a, b, x, 0.0, normalised, ref p_derivative, y);
                            }
                            else
                            {
                                fract = -(normalised ? 1 : beta(a, b));
                                invert = false;
                                fract = -ibeta_series(a, b, x, fract, normalised, ref p_derivative, y);
                            }
                        }
                        else if (a >= 15)
                        {
                            if (!invert)
                            {
                                fract = beta_small_b_large_a_series(a, b, x, y, 0.0, 1.0, normalised);
                            }
                            else
                            {
                                fract = -(normalised ? 1 : beta(a, b));
                                invert = false;
                                fract = -beta_small_b_large_a_series(a, b, x, y, fract, 1.0, normalised);
                            }
                        }
                        else
                        {
                            // Sidestep to improve errors:
                            double prefix;
                            if (!normalised)
                            {
                                prefix = rising_factorial_ratio(a + b, a, 20);
                            }
                            else
                            {
                                prefix = 1;
                            }
                            fract = ibeta_a_step(a, b, x, y, 20, normalised, ref p_derivative);
                            if (!invert)
                            {
                                fract = beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, normalised);
                            }
                            else
                            {
                                fract -= (normalised ? 1 : beta(a, b));
                                invert = false;
                                fract = -beta_small_b_large_a_series(a + 20, b, x, y, fract, prefix, normalised);
                            }
                        }
                    }
                }
            }
            else
            {
                // Both a,b >= 1:
                double lambda;
                if (a < b)
                {
                    lambda = a - (a + b) * x;
                }
                else
                {
                    lambda = (a + b) * y - b;
                }
                if (lambda < 0)
                {
                    swap(ref a, ref b);
                    swap(ref x, ref y);
                    invert = !invert;
                }

                if (b < 40)
                {
                    if ((Math.Floor(a) == a) && (Math.Floor(b) == b))
                    {
                        // relate to the binomial distribution and use a finite sum:
                        double k = a - 1;
                        double n = b + k;
                        fract = binomial_ccdf(n, k, x, y);
                        if (!normalised)
                            fract *= beta(a, b);
                    }
                    else if (b * x <= 0.7)
                    {
                        if (!invert)
                        {
                            fract = ibeta_series(a, b, x, 0.0, normalised, ref p_derivative, y);
                        }
                        else
                        {
                            fract = -(normalised ? 1 : beta(a, b));
                            invert = false;
                            fract = -ibeta_series(a, b, x, fract, normalised, ref p_derivative, y);
                        }
                    }
                    else if (a > 15)
                    {
                        // sidestep so we can use the series representation:
                        int n = (int)Math.Floor(b);
                        if (n == b)
                            --n;
                        double bbar = b - n;
                        double prefix;
                        if (!normalised)
                        {
                            prefix = rising_factorial_ratio(a + bbar, bbar, n);
                        }
                        else
                        {
                            prefix = 1;
                        }
                        double pd = 0.0;
                        fract = ibeta_a_step(bbar, a, y, x, n, normalised, ref pd);
                        fract = beta_small_b_large_a_series(a, bbar, x, y, fract, 1.0, normalised);
                        fract /= prefix;
                    }
                    else if (normalised)
                    {
                        // the formula here for the non-normalised case is tricky to figure
                        // out (for me!!), and requires two pochhammer calculations rather
                        // than one, so leave it for now....
                        int n = (int)Math.Floor(b);
                        double bbar = b - n;
                        if (bbar <= 0)
                        {
                            --n;
                            bbar += 1;
                        }
                        double pd = 0.0;
                        fract = ibeta_a_step(bbar, a, y, x, n, normalised, ref pd);
                        fract += ibeta_a_step(a, bbar, x, y, 20, normalised, ref pd);
                        if (invert)
                            fract -= (normalised ? 1 : beta(a, b));
                        //fract = ibeta_series(a+20, bbar, x, fract, l, normalised, p_derivative, y);
                        fract = beta_small_b_large_a_series(a + 20, bbar, x, y, fract, 1.0, normalised);
                        if (invert)
                        {
                            fract = -fract;
                            invert = false;
                        }
                    }
                    else
                    {
                        fract = ibeta_fraction2(a, b, x, y, normalised, ref p_derivative);
                    }
                }
                else
                {
                    fract = ibeta_fraction2(a, b, x, y, normalised, ref p_derivative);
                }
            }

            if (p_derivative < 0) p_derivative = ibeta_power_terms(a, b, x, y, true);
            double div = y * x;

            if (p_derivative != 0)
            {
                if ((double.MaxValue * div < p_derivative))
                {
                    // overflow, return an arbitarily large value:
                    p_derivative = double.MaxValue / 2;
                }
                else
                {
                    p_derivative /= div;
                }
            }
            return invert ? (normalised ? 1 : beta(a, b)) - fract : fract;
        }

        private static double ibeta_fraction2(double a, double b, double x, double y, bool normalised, ref double p_derivative)
        {
            double result = ibeta_power_terms(a, b, x, y, normalised);
            p_derivative = result;
            if (result == 0) return result;

            ibeta_fraction2_t f = new ibeta_fraction2_t(a, b, x);
            int terms = max_series_iterations;
            double fract = continued_fraction_b(f, XMath.epsilon, ref terms);
            return result / fract;
        }

        class ibeta_fraction2_t : series<pair<double>>
        {
            double a, b, x;
            int m;

            public ibeta_fraction2_t(double a_, double b_, double x_)
            {
                a = a_;
                b = b_;
                x = x_;
                m = 0;
            }

            public override pair<double> next()
            {
                double aN = (a + m - 1) * (a + b + m - 1) * m * (b - m) * x * x;
                double denom = (a + 2 * m - 1);
                aN /= denom * denom;

                double bN = m;
                bN += (m * (b - m) * x) / (a + 2 * m - 1);
                bN += ((a + m) * (a - (a + b) * x + 1 + m * (2 - x))) / (a + 2 * m + 1);

                ++m;

                pair<double> p = new pair<double>();
                p.v1 = aN;
                p.v2 = bN;
                return p;
            }
        }

        internal static double ibeta_inv_imp(double a, double b, double p, double q, ref double py)
        {
            //
            // The flag invert is set to true if we swap a for b and p for q,
            // in which case the result has to be subtracted from 1:
            //
            bool invert = false;
            //
            // Depending upon which approximation method we use, we may end up
            // calculating either x or y initially (where y = 1-x):
            //
            double x = 0; // Set to a safe zero to avoid a
            // MSVC 2005 warning C4701: potentially uninitialized local variable 'x' used
            // But code inspection appears to ensure that x IS assigned whatever the code path.
            double y = 0.0;

            // For some of the methods we can put tighter bounds
            // on the result than simply [0,1]:
            //
            double lower = 0;
            double upper = 1;
            //
            // Student's double with b = 0.5 gets handled as a special case, swap
            // around if the arguments are in the "wrong" order:
            //
            if (a == 0.5)
            {
                swap(ref a, ref b);
                swap(ref p, ref q);
                invert = !invert;
            }
            //
            // Handle trivial cases first:
            //
            if (q == 0)
            {
                py = 0;
                return 1;
            }
            else if (p == 0)
            {
                py = 1;
                return 0;
            }
            else if ((a == 1) && (b == 1))
            {
                py = 1 - p;
                return p;
            }
            else if ((b == 0.5) && (a >= 0.5))
            {
                //
                // We have a Student's double distribution:
                x = Distributions.Students_T_Distribution.find_ibeta_inv_from_t_dist(a, p, q, ref y);
            }
            else if (a + b > 5)
            {
                //
                // When a+b is large then we can use one of Prof Temme's
                // asymptotic expansions, begin by swapping things around
                // so that p < 0.5, we do this to avoid cancellations errors
                // when p is large.
                //
                if (p > 0.5)
                {
                    swap(ref a, ref b);
                    swap(ref p, ref q);
                    invert = !invert;
                }
                double minv = min(a, b);
                double maxv = max(a, b);
                if ((Math.Sqrt(minv) > (maxv - minv)) && (minv > 5))
                {
                    //
                    // When a and b differ by a small amount
                    // the curve is quite symmetrical and we can use an error
                    // function to approximate the inverse. This is the cheapest
                    // of the three Temme expantions, and the calculated value
                    // for x will never be much larger than p, so we don't have
                    // to worry about cancellation as long as p is small.
                    //
                    x = temme_method_1_ibeta_inverse(a, b, p);
                    y = 1 - x;
                }
                else
                {
                    double r = a + b;
                    double theta = Math.Asin(Math.Sqrt(a / r));
                    double lambda = minv / r;
//                    if ((lambda >= 0.2) && (lambda <= 0.8) && (lambda >= 10))     //CGS Change - this line could never evaluate to true !! Must be a typo !
                    if ((lambda >= 0.2) && (lambda <= 0.8) && (r >= 10))       //I'm guessing here...
                    {
                        //
                        // The second error function case is the next cheapest
                        // to use, it brakes down when the result is likely to be
                        // very small, if a+b is also small, but we can use a
                        // cheaper expansion there in any case.  As before x won't
                        // be much larger than p, so as long as p is small we should
                        // be free of cancellation error.
                        //
                        double ppa = Math.Pow(p, 1 / a);
                        if ((ppa < 0.0025) && (a + b < 200))
                        {
                            x = ppa * Math.Pow(a * beta(a, b), 1 / a);
                        }
                        else
                            x = temme_method_2_ibeta_inverse(a, b, p, r, theta);
                        y = 1 - x;
                    }
                    else
                    {
                        //
                        // If we get here then a and b are very different in magnitude
                        // and we need to use the third of Temme's methods which
                        // involves inverting the incomplete gamma.  This is much more
                        // expensive than the other methods.  We also can only use this
                        // method when a > b, which can lead to cancellation errors
                        // if we really want y (as we will when x is close to 1), so
                        // a different expansion is used in that case.
                        //
                        if (a < b)
                        {
                            swap(ref a, ref b);
                            swap(ref p, ref q);
                            invert = !invert;
                        }
                        //
                        // Try and compute the easy way first:
                        //
                        double bet = 0;
                        if (b < 2)
                            bet = beta(a, b);
                        if (bet != 0)
                        {
                            y = Math.Pow(b * q * bet, 1 / b);
                            x = 1 - y;
                        }
                        else
                            y = 1;
                        if (y > 1e-5)
                        {
                            x = temme_method_3_ibeta_inverse(a, b, p, q);
                            y = 1 - x;
                        }
                    }
                }
            }
            else if ((a < 1) && (b < 1))
            {
                //
                // Both a and b less than 1,
                // there is a point of inflection at xs:
                //
                double xs = (1 - a) / (2 - a - b);
                //
                // Now we need to ensure that we start our iteration from the
                // right side of the inflection point:
                //
                double fs = ibeta(a, b, xs) - p;
                if (Math.Abs(fs) / p < XMath.epsilon * 3)
                {
                    // The result is at the point of inflection, best just return it:
                    py = invert ? xs : 1 - xs;
                    return invert ? 1 - xs : xs;
                }
                if (fs < 0)
                {
                    swap(ref a, ref b);
                    swap(ref p, ref q);
// CGS CHANGE - ATTEMPT TO FIX BUG
//                    invert = true;
                    invert = !invert;       
// END
                    xs = 1 - xs;
                }
                double xg = Math.Pow(a * p * beta(a, b), 1 / a);
                x = xg / (1 + xg);
                y = 1 / (1 + xg);
                //
                // And finally we know that our result is below the inflection
                // point, so set an upper limit on our search:
                //
                if (x > xs)
                    x = xs;
                upper = xs;
            }
            else if ((a > 1) && (b > 1))
            {
                //
                // Small a and b, both greater than 1,
                // there is a point of inflection at xs,
                // and it's complement is xs2, we must always
                // start our iteration from the right side of the
                // point of inflection.
                //
                double xs = (a - 1) / (a + b - 2);
                double xs2 = (b - 1) / (a + b - 2);
                double ps = ibeta(a, b, xs) - p;

                if (ps < 0)
                {
                    swap(ref a, ref b);
                    swap(ref p, ref q);
                    swap(ref xs, ref xs2);
                    invert = true;
                }
                //
                // Estimate x and y, using expm1 to get a good estimate
                // for y when it's very small:
                //
                double lx = Math.Log(p * a * beta(a, b)) / a;
                x = Math.Exp(lx);
                y = x < 0.9 ? 1.0 - x : (double)(-expm1(lx));

                if ((b < a) && (x < 0.2))
                {
                    //
                    // Under a limited range of circumstances we can improve
                    // our estimate for x, frankly it's not clear if this has much effect!
                    //
                    double ap1 = a - 1;
                    double bm1 = b - 1;
                    double a_2 = a * a;
                    double a_3 = a * a_2;
                    double b_2 = b * b;
                    double[] terms = new double[5];
                    terms[0] = 0;
                    terms[1] = 1;
                    terms[2] = bm1 / ap1;
                    ap1 *= ap1;
                    terms[3] = bm1 * (3 * a * b + 5 * b + a_2 - a - 4) / (2 * (a + 2) * ap1);
                    ap1 *= (a + 1);
                    terms[4] = bm1 * (33 * a * b_2 + 31 * b_2 + 8 * a_2 * b_2 - 30 * a * b - 47 * b + 11 * a_2 * b + 6 * a_3 * b + 18 + 4 * a - a_3 + a_2 * a_2 - 10 * a_2)
                               / (3 * (a + 3) * (a + 2) * ap1);
                    x = evaluate_polynomial(terms, x, 5);
                }
                //
                // And finally we know that our result is below the inflection
                // point, so set an upper limit on our search:
                //
                if (x > xs)
                    x = xs;
                upper = xs;
            }
            else /*if((a <= 1) != (b <= 1))*/
            {
                //
                // If all else fails we get here, only one of a and b
                // is above 1, and a+b is small.  Start by swapping
                // things around so that we have a concave curve with b > a
                // and no points of inflection in [0,1].  As long as we expect
                // x to be small then we can use the simple (and cheap) power
                // term to estimate x, but when we expect x to be large then
                // this greatly underestimates x and leaves us trying to
                // iterate "round the corner" which may take almost forever...
                //
                // We could use Temme's inverse gamma function case in that case,
                // this works really rather well (albeit expensively) even though
                // strictly speaking we're outside it's defined range.
                //
                // However it's expensive to compute, and an alternative approach
                // which models the curve as a distorted quarter circle is much
                // cheaper to compute, and still keeps the number of iterations
                // required down to a reasonable level.  With thanks to Prof Temme
                // for this suggestion.
                //
                if (b < a)
                {
                    swap(ref a, ref b);
                    swap(ref p, ref q);
                    invert = true;
                }
                if (Math.Pow(p, 1 / a) < 0.5)
                {
                    x = Math.Pow(p * a * beta(a, b), 1 / a);
                    if (x == 0)
                        x = XMath.min_value;
                    y = 1 - x;
                }
                else /*if(Math.Pow(q, 1/b) < 0.1)*/
                {
                    // model a distorted quarter circle:
                    y = Math.Pow(1 - Math.Pow(p, b * beta(a, b)), 1 / b);
                    if (y == 0)
                        y = XMath.min_value;
                    x = 1 - y;
                }
            }

            //
            // Now we have a guess for x (and for y) we can set things up for
            // iteration.  If x > 0.5 it pays to swap things round:
            //
            if (x > 0.5)
            {
                swap(ref a, ref b);
                swap(ref p, ref q);
                swap(ref x, ref y);
                invert = !invert;
                double l = 1 - upper;
                double u = 1 - lower;
                lower = l;
                upper = u;
            }
            //
            // lower bound for our search:
            //
            // We're not interested in denormalised answers as these tend to
            // these tend to take up lots of iterations, given that we can't get
            // accurate derivatives in this area (they tend to be infinite).
            //
            if (lower == 0)
            {
                if (invert && (py == 0))
                {
                    //
                    // We're not interested in answers smaller than machine epsilon:
                    //
                    lower = double.Epsilon;
                    if (x < lower)
                        x = lower;
                }
                else
                    lower = double.Epsilon;
                if (x < lower)
                    x = lower;
            }
            //
            // Figure out how many digits to iterate towards:
            //
            int digits = 53/2;
            if ((x < 1e-50) && ((a < 1) || (b < 1)))
            {
                //
                // If we're in a region where the first derivative is very
                // large, then we have to take care that the root-finder
                // doesn't terminate prematurely.  We'll bump the precision
                // up to avoid this, but we have to take care not to set the
                // precision too high or the last few iterations will just
                // thrash around and convergence may be slow in this case.
                // Try 3/4 of machine epsilon:
                //
                digits *= 3;
                digits /= 2;
            }
            //
            // Now iterate, we can use either p or q as the target here
            // depending on which is smaller:
            //
            int max_iter = max_root_iterations;
            x = halley_iterate(new ibeta_roots(a, b, (p < q ? p : q), (p < q ? false : true)), x, lower, upper, digits, max_iter);
            //
            // Tidy up, if "lower" was too high then zero is the best answer we have:
            //
            if (x == lower) x = 0;
            py = invert ? x : 1 - x;
            return invert ? 1 - x : x;
        }

        private class beta_inv_ab_t : iterand<double, double>
        {
            double b, z, p;
            bool invert, swap_ab;

            public beta_inv_ab_t(double b_, double z_, double p_, bool invert_, bool swap_ab_)
            {
                b = b_;
                z = z_;
                p = p_;
                invert = invert_;
                swap_ab = swap_ab_;
            }
            public override double next(double a)
            {
                return invert ?
                   p - ibetac(swap_ab ? b : a, swap_ab ? a : b, z) : ibeta(swap_ab ? b : a, swap_ab ? a : b, z) - p;
            }
        }

        private static double ibeta_inv_ab_imp(double b, double z, double p, double q, bool swap_ab)
        {
            if (p == 0) return swap_ab ? XMath.min_value : double.MaxValue;
            if (q == 0) return swap_ab ? double.MaxValue : XMath.min_value;
            //
            // Function object, this is the functor whose root
            // we have to solve:
            //
            beta_inv_ab_t f = new beta_inv_ab_t(b, z, (p < q) ? p : q, (p < q) ? false : true, swap_ab);
            //
            // Tolerance: full precision.
            //
            eps_tolerance tol = new eps_tolerance(53);
            //
            // Now figure out a starting guess for what a may be, 
            // we'll start out with a value that'll put p or q
            // right bang in the middle of their range, the functions
            // are quite sensitive so we should need too many steps
            // to bracket the root from there:
            //
            double guess = 0;
            double factor = 5;
            //
            // Convert variables to parameters of a negative binomial distribution:
            //
            double n = b;
            double sf = swap_ab ? z : 1 - z;
            double sfc = swap_ab ? 1 - z : z;
            double u = swap_ab ? p : q;
            double v = swap_ab ? q : p;
            if (u <= Math.Pow(sf, n))
            {
                //
                // Result is less than 1, negative binomial approximation
                // is useless....
                //
                if ((p < q) != swap_ab) guess = min(b * 2.0, 1.0);
                else guess = min(b / 2.0, 1.0);
            }
            if (n * n * n * u * sf > 0.005) guess = 1 + inverse_negative_binomial_cornish_fisher(n, sf, sfc, u, v);

            if (guess < 10)
            {
                //
                // Negative binomial approximation not accurate in this area:
                //
                if ((p < q) != swap_ab) guess = min(b * 2.0, 10.0);
                else guess = min(b / 2.0, 10.0);
            }
            else
                factor = (v < Math.Sqrt(XMath.epsilon)) ? 2 : (guess < 20 ? 1.2 : 1.1);
            //
            // Max iterations permitted:
            //
            int max_iter = max_root_iterations;
            pair<double> r = bracket_and_solve_root(f, guess, factor, swap_ab ? true : false, tol, ref max_iter);
            if (max_iter >= max_root_iterations)
                throw new Exception(string.Format("ibeta_inv_ab: unable to locate the root within a reasonable number of iterations, closest approximation so far was {0:G}", r.v1));
            return (r.v1 + r.v2) / 2;
        }

        //
        // See:
        // "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157.
        // Section 2.
        //
        private static double temme_method_1_ibeta_inverse(double a, double b, double z)
        {
            double r2 = XMath.root_two;
            //
            // get the first approximation for eta from the inverse
            // error function (Eq: 2.9 and 2.10).
            //
            double eta0 = erfc_inv(2 * z);
            eta0 /= -Math.Sqrt(a / 2);

            double[] terms = new double[4];
            terms[0] = eta0;
            double[] workspace = new double[7];
            //
            // calculate powers:
            //
            double B = b - a;
            double B_2 = B * B;
            double B_3 = B_2 * B;
            //
            // Calculate correction terms:
            //

            // See eq following 2.15:
            workspace[0] = -B * r2 / 2;
            workspace[1] = (1 - 2 * B) / 8;
            workspace[2] = -(B * r2 / 48);
            workspace[3] = -1.0 / 192.0;
            workspace[4] = -B * r2 / 3840;
            terms[1] = evaluate_polynomial(workspace, eta0, 5);
            // Eq Following 2.17:
            workspace[0] = B * r2 * (3 * B - 2) / 12;
            workspace[1] = (20 * B_2 - 12 * B + 1) / 128;
            workspace[2] = B * r2 * (20 * B - 1) / 960;
            workspace[3] = (16 * B_2 + 30 * B - 15) / 4608;
            workspace[4] = B * r2 * (21 * B + 32) / 53760;
            workspace[5] = (-32 * B_2 + 63) / 368640;
            workspace[6] = -B * r2 * (120 * B + 17) / 25804480;
            terms[2] = evaluate_polynomial(workspace, eta0, 7);
            // Eq Following 2.17:
            workspace[0] = B * r2 * (-75 * B_2 + 80 * B - 16) / 480;
            workspace[1] = (-1080 * B_3 + 868 * B_2 - 90 * B - 45) / 9216;
            workspace[2] = B * r2 * (-1190 * B_2 + 84 * B + 373) / 53760;
            workspace[3] = (-2240 * B_3 - 2508 * B_2 + 2100 * B - 165) / 368640;
            terms[3] = evaluate_polynomial(workspace, eta0, 4);
            //
            // Bring them together to get a final estimate for eta:
            //
            double eta = evaluate_polynomial(terms, 1.0 / a, 4);
            //
            // now we need to convert eta to x, by solving the appropriate
            // quadratic equation:
            //
            double eta_2 = eta * eta;
            double c = -Math.Exp(-eta_2 / 2);
            double x;
            if (eta_2 == 0)
                x = 0.5;
            else
                x = (1 + eta * Math.Sqrt((1 + c) / eta_2)) / 2;

            return x;
        }
        //
        // See:
        // "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157.
        // Section 3.
        //
        private static double temme_method_2_ibeta_inverse(double a, double b, double z, double r, double theta)
        {
            //
            // Get first estimate for eta, see Eq 3.9 and 3.10,
            // but note there is a typo in Eq 3.10:
            //
            double eta0 = erfc_inv(2 * z);
            eta0 /= -Math.Sqrt(r / 2);

            double s = Math.Sin(theta);
            double c = Math.Cos(theta);
            //
            // Now we need to purturb eta0 to get eta, which we do by
            // evaluating the polynomial in 1/r at the bottom of page 151,
            // to do this we first need the error terms e1, e2 e3
            // which we'll fill into the array "terms".  Since these
            // terms are themselves polynomials, we'll need another
            // array "workspace" to calculate those...
            //
            double[] terms = new double[4];
            terms[0] = eta0;
            double[] workspace = new double[6];
            //
            // some powers of sin(theta)cos(theta) that we'll need later:
            //
            double sc = s * c;
            double sc_2 = sc * sc;
            double sc_3 = sc_2 * sc;
            double sc_4 = sc_2 * sc_2;
            double sc_5 = sc_2 * sc_3;
            double sc_6 = sc_3 * sc_3;
            double sc_7 = sc_4 * sc_3;
            //
            // Calculate e1 and put it in terms[1], see the middle of page 151:
            //
            workspace[0] = (2 * s * s - 1) / (3 * s * c);
            double[] co1 = { -1, -5, 5 };
            workspace[1] = -evaluate_even_polynomial(co1, s, 3) / (36 * sc_2);
            double[] co2 = { 1, 21, -69, 46 };
            workspace[2] = evaluate_even_polynomial(co2, s, 4) / (1620 * sc_3);
            double[] co3 = { 7, -2, 33, -62, 31 };
            workspace[3] = -evaluate_even_polynomial(co3, s, 5) / (6480 * sc_4);
            double[] co4 = { 25, -52, -17, 88, -115, 46 };
            workspace[4] = evaluate_even_polynomial(co4, s, 6) / (90720 * sc_5);
            terms[1] = evaluate_polynomial(workspace, eta0, 5);
            //
            // Now evaluate e2 and put it in terms[2]:
            //
            double[] co5 = { 7, 12, -78, 52 };
            workspace[0] = -evaluate_even_polynomial(co5, s, 4) / (405 * sc_3);
            double[] co6 = { -7, 2, 183, -370, 185 };
            workspace[1] = evaluate_even_polynomial(co6, s, 5) / (2592 * sc_4);
            double[] co7 = { -533, 776, -1835, 10240, -13525, 5410 };
            workspace[2] = -evaluate_even_polynomial(co7, s, 6) / (204120 * sc_5);
            double[] co8 = { -1579, 3747, -3372, -15821, 45588, -45213, 15071 };
            workspace[3] = -evaluate_even_polynomial(co8, s, 7) / (2099520 * sc_6);
            terms[2] = evaluate_polynomial(workspace, eta0, 4);
            //
            // And e3, and put it in terms[3]:
            //
            double[] co9 = {449, -1259, -769, 6686, -9260, 3704 };
            workspace[0] = evaluate_even_polynomial(co9, s, 6) / (102060 * sc_5);
            double[] co10 = { 63149, -151557, 140052, -727469, 2239932, -2251437, 750479 };
            workspace[1] = -evaluate_even_polynomial(co10, s, 7) / (20995200 * sc_6);
            double[] co11 = { 29233, -78755, 105222, 146879, -1602610, 3195183, -2554139, 729754 };
            workspace[2] = evaluate_even_polynomial(co11, s, 8) / (36741600 * sc_7);
            terms[3] = evaluate_polynomial(workspace, eta0, 3);
            //
            // Bring the correction terms together to evaluate eta,
            // this is the last equation on page 151:
            //
            double eta = evaluate_polynomial(terms, 1.0/r, 4);
            //
            // Now that we have eta we need to back solve for x,
            // we seek the value of x that gives eta in Eq 3.2.
            // The two methods used are described in section 5.
            //
            // Begin by defining a few variables we'll need later:
            //
            double x;
            double s_2 = s * s;
            double c_2 = c * c;
            double alpha = c / s;
            alpha *= alpha;
            double lu = (-(eta * eta) / (2 * s_2) + Math.Log(s_2) + c_2 * Math.Log(c_2) / s_2);
            //
            // Temme doesn't specify what value to switch on here,
            // but this seems to work pretty well:
            //
            if(Math.Abs(eta) < 0.7)
            {
                //
                // Small eta use the expansion Temme gives in the second equation
                // of section 5, it's a polynomial in eta:
                //
                workspace[0] = s * s;
                workspace[1] = s * c;
                workspace[2] = (1 - 2 * workspace[0]) / 3;
                double[] co12 = { 1, -13, 13 };
                workspace[3] = evaluate_polynomial(co12, workspace[0], 3) / (36 * s * c);
                double[] co13 = { 1, 21, -69, 46 };
                workspace[4] = evaluate_polynomial(co13, workspace[0], 4) / (270 * workspace[0] * c * c);
                x = evaluate_polynomial(workspace, eta, 5);
            }
            else
            {
                //
                // If eta is large we need to solve Eq 3.2 more directly,
                // begin by getting an initial approximation for x from
                // the last equation on page 155, this is a polynomial in u:
                //
                double u = Math.Exp(lu);
                workspace[0] = u;
                workspace[1] = alpha;
                workspace[2] = 0;
                workspace[3] = 3 * alpha * (3 * alpha + 1) / 6;
                workspace[4] = 4 * alpha * (4 * alpha + 1) * (4 * alpha + 2) / 24;
                workspace[5] = 5 * alpha * (5 * alpha + 1) * (5 * alpha + 2) * (5 * alpha + 3) / 120;
                x = evaluate_polynomial(workspace, u, 6);
                //
                // At this point we may or may not have the right answer, Eq-3.2 has
                // two solutions for x for any given eta, however the mapping in 3.2
                // is 1:1 with the sign of eta and x-sin^2(theta) being the same.
                // So we can check if we have the right root of 3.2, and if not
                // switch x for 1-x.  This transformation is motivated by the fact
                // that the distribution is *almost* symetric so 1-x will be in the right
                // ball park for the solution:
                //
                if((x - s_2) * eta < 0)
                 x = 1 - x;
            }
            //
            // The final step is a few Newton-Raphson iterations to
            // clean up our approximation for x, this is pretty cheap
            // in general, and very cheap compared to an incomplete beta
            // evaluation.  The limits set on x come from the observation
            // that the sign of eta and x-sin^2(theta) are the same.
            //
            double lower, upper;
            if(eta < 0)
            {
                lower = 0;
                upper = s_2;
            }
            else
            {
                lower = s_2;
                upper = 1;
            }
            //
            // If our initial approximation is out of bounds then bisect:
            //
            if((x < lower) || (x > upper)) x = (lower+upper) / 2;
            //
            // And iterate:
            //
            x = newton_raphson_iterate(new temme_root_finder(-lu, alpha), x, lower, upper, 53 / 2, max_root_iterations);

            return x;
        }
        //
        // See:
        // "Asymptotic Inversion of the Incomplete Beta Function"
        // N.M. Temme
        // Journal of Computation and Applied Mathematics 41 (1992) 145-157.
        // Section 4.
        //
        private static double temme_method_3_ibeta_inverse(double a, double b, double p, double q)
        {
            //
            // Begin by getting an initial approximation for the quantity
            // eta from the dominant part of the incomplete beta:
            //
            double eta0;
            if (p < q)
                eta0 = gamma_q_inv(b, p);
            else
                eta0 = gamma_p_inv(b, q);
            eta0 /= a;
            //
            // Define the variables and powers we'll need later on:
            //
            double mu = b / a;
            double w = Math.Sqrt(1 + mu);
            double w_2 = w * w;
            double w_3 = w_2 * w;
            double w_4 = w_2 * w_2;
            double w_5 = w_3 * w_2;
            double w_6 = w_3 * w_3;
            double w_7 = w_4 * w_3;
            double w_8 = w_4 * w_4;
            double w_9 = w_5 * w_4;
            double w_10 = w_5 * w_5;
            double d = eta0 - mu;
            double d_2 = d * d;
            double d_3 = d_2 * d;
            double d_4 = d_2 * d_2;
            double w1 = w + 1;
            double w1_2 = w1 * w1;
            double w1_3 = w1 * w1_2;
            double w1_4 = w1_2 * w1_2;
            //
            // Now we need to compute the purturbation error terms that
            // convert eta0 to eta, these are all polynomials of polynomials.
            // Probably these should be re-written to use tabulated data
            // (see examples above), but it's less of a win in this case as we
            // need to calculate the individual powers for the denominator terms
            // anyway, so we might as well use them for the numerator-polynomials
            // as well....
            //
            // Refer to p154-p155 for the details of these expansions:
            //
            double e1 = (w + 2) * (w - 1) / (3 * w);
            e1 += (w_3 + 9 * w_2 + 21 * w + 5) * d / (36 * w_2 * w1);
            e1 -= (w_4 - 13 * w_3 + 69 * w_2 + 167 * w + 46) * d_2 / (1620 * w1_2 * w_3);
            e1 -= (7 * w_5 + 21 * w_4 + 70 * w_3 + 26 * w_2 - 93 * w - 31) * d_3 / (6480 * w1_3 * w_4);
            e1 -= (75 * w_6 + 202 * w_5 + 188 * w_4 - 888 * w_3 - 1345 * w_2 + 118 * w + 138) * d_4 / (272160 * w1_4 * w_5);

            double e2 = (28 * w_4 + 131 * w_3 + 402 * w_2 + 581 * w + 208) * (w - 1) / (1620 * w1 * w_3);
            e2 -= (35 * w_6 - 154 * w_5 - 623 * w_4 - 1636 * w_3 - 3983 * w_2 - 3514 * w - 925) * d / (12960 * w1_2 * w_4);
            e2 -= (2132 * w_7 + 7915 * w_6 + 16821 * w_5 + 35066 * w_4 + 87490 * w_3 + 141183 * w_2 + 95993 * w + 21640) * d_2 / (816480 * w_5 * w1_3);
            e2 -= (11053 * w_8 + 53308 * w_7 + 117010 * w_6 + 163924 * w_5 + 116188 * w_4 - 258428 * w_3 - 677042 * w_2 - 481940 * w - 105497) * d_3 / (14696640 * w1_4 * w_6);

            double e3 = -((3592 * w_7 + 8375 * w_6 - 1323 * w_5 - 29198 * w_4 - 89578 * w_3 - 154413 * w_2 - 116063 * w - 29632) * (w - 1)) / (816480 * w_5 * w1_2);
            e3 -= (442043 * w_9 + 2054169 * w_8 + 3803094 * w_7 + 3470754 * w_6 + 2141568 * w_5 - 2393568 * w_4 - 19904934 * w_3 - 34714674 * w_2 - 23128299 * w - 5253353) * d / (146966400 * w_6 * w1_3);
            e3 -= (116932 * w_10 + 819281 * w_9 + 2378172 * w_8 + 4341330 * w_7 + 6806004 * w_6 + 10622748 * w_5 + 18739500 * w_4 + 30651894 * w_3 + 30869976 * w_2 + 15431867 * w + 2919016) * d_2 / (146966400 * w1_4 * w_7);
            //
            // Combine eta0 and the error terms to compute eta (Second eqaution p155):
            //
            double eta = eta0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
            //
            // Now we need to solve Eq 4.2 to obtain x.  For any given value of
            // eta there are two solutions to this equation, and since the distribtion
            // may be very skewed, these are not related by x ~ 1-x we used when
            // implementing section 3 above.  However we know that:
            //
            //  cross < x <= 1       ; iff eta < mu
            //          x == cross   ; iff eta == mu
            //     0 <= x < cross    ; iff eta > mu
            //
            // Where cross == 1 / (1 + mu)
            // Many thanks to Prof Temme for clarifying this point.
            //
            // Therefore we'll just jump straight into Newton iterations
            // to solve Eq 4.2 using these bounds, and simple bisection
            // as the first guess, in practice this converges pretty quickly
            // and we only need a few digits correct anyway:
            //
            if (eta <= 0)
                eta = XMath.min_value;
            double u = eta - mu * Math.Log(eta) + (1 + mu) * Math.Log(1 + mu) - mu;
            double cross = 1 / (1 + mu);
            double lower = eta < mu ? cross : 0;
            double upper = eta < mu ? 1 : cross;
            double x = (lower + upper) / 2;
            x = newton_raphson_iterate(new temme_root_finder(u, mu), x, lower, upper, 53 / 2, max_root_iterations);
            return x;
        }

        class temme_root_finder : iterand<Tuple<double, double>, double>
        {
            double t, a;

            public temme_root_finder(double t_, double a_)
            {
                t = t_;
                a = a_;
            }

            public override Tuple<double, double> next(double x)
            {
                double y = 1 - x;
                if (y == 0)
                {
                    double big = double.MaxValue / 4;
                    return new Tuple<double, double>(-big, -big);
                }
                if (x == 0)
                {
                    double big = double.MaxValue / 4;
                    return new Tuple<double, double>(-big, big);
                }
                double f = Math.Log(x) + a * Math.Log(y) + t;
                double f1 = (1 / x) - (a / (y));
                return new Tuple<double, double>(f, f1);
            }
        }

        class ibeta_roots : iterand<Tuple<double, double, double>, double>
        {
            double a, b, target;
            bool invert;

            public ibeta_roots(double _a, double _b, double t, bool inv)
            {
                a = _a;
                b = _b;
                target = t;
                invert = inv;
            }

            public override Tuple<double, double, double> next(double x)
            {
                double f1 = 0.0;
                double y = 1 - x;
                double f = ibeta_imp(a, b, x, invert, true, ref f1) - target;
                if (invert)
                    f1 = -f1;
                if (y == 0)
                    y = XMath.min_value * 64;
                if (x == 0)
                    x = XMath.min_value * 64;

                double f2 = f1 * (-y * a + (b - 2) * x + 1);
                if (Math.Abs(f2) < y * x * double.MaxValue)
                    f2 /= (y * x);
                if (invert)
                    f2 = -f2;

                // make sure we don't have a zero derivative:
                if (f1 == 0)
                    f1 = (invert ? -1 : 1) * double.Epsilon * 64;

                return new Tuple<double, double, double>(f, f1, f2);
            }
        }

        #endregion
    }
}
