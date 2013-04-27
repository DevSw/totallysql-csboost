using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost
{
    public static partial class XMath
    {
        public static double gamma(double z)
        {
            double result = 1;

            if (z <= 0)
            {
                if (Math.Floor(z) == z) throw new Exception("Evaluation of tgamma at a negative integer " + z.ToString());
                if (z <= -20.0)
                {
                    result = gamma(-z) * sinpx(z);
                    if ((Math.Abs(result) < 1) && (double.MaxValue * Math.Abs(result) < Math.PI)) throw new Exception("Result of tgamma is too large to represent.");
                    result = -Math.PI / result;
                    if (result == 0) throw new Exception("Result of gamma is too small to represent.");
                    return result;
                }
                // shift z to > 1:
                while (z < 0)
                {
                    result /= z;
                    z += 1;
                }
            }
            double[] f = factorials();
            if ((Math.Floor(z) == z) && (z < f.Length))
            {
                result *= f[(int)z - 1];
            }
            else
            {
                result *= lanczos_sum(z);
                if (z * Math.Log(z) > log_max_value)
                {
                    // we're going to overflow unless this is done with care:
                    double zgh = (z + lanczos_g - 0.5);
                    if (Math.Log(zgh) * z / 2 > log_max_value) throw new Exception("Result of tgamma is too large to represent.");
                    double hp = Math.Pow(zgh, (z / 2.0) - 0.25);
                    result *= hp / Math.Exp(zgh);
                    if (double.MaxValue / hp < result) throw new Exception("Result of tgamma is too large to represent.");
                    result *= hp;
                }
                else
                {
                    double zgh = (z + lanczos_g - 0.5);
                    result *= Math.Pow(zgh, z - 0.5) / Math.Exp(zgh);
                }
            }
            return result;
        }

        public static double gamma_p(double a, double z)
        {
            double p_derivative = 0.0;
            return gamma_incomplete_imp(a, z, true, false, ref p_derivative);
        }

        public static double gamma_lower(double a, double z)
        {
            double p_derivative = 0.0;
            return gamma_incomplete_imp(a, z, false, false, ref p_derivative);
        }

        public static double gamma_q(double a, double x)
        {
            double p_derivative = 0.0;
            return gamma_incomplete_imp(a, x, true, true, ref p_derivative);
        }

        public static double gamma_upper(double a, double x)
        {
            double p_derivative = 0.0;
            return gamma_incomplete_imp(a, x, false, true, ref p_derivative);
        }

        public static double gamma_q_inv(double a, double q)
        {
            if (a <= 0)
                throw new Exception(string.Format("Argument a in the incomplete gamma function inverse must be >= 0 (got a={0:G}).", a));
            if ((q < 0) || (q > 1))
                throw new Exception(string.Format("Probabilty must be in the range [0,1] in the incomplete gamma function inverse (got q={0:G}).", q));
            if (q == 0)
                return double.MaxValue;
            if (q == 1)
                return 0;
            bool has_10_digits;
            double guess = find_inverse_gamma(a, 1 - q, q, out has_10_digits);
            double lower = XMath.min_value;
            if (guess <= lower) guess = XMath.min_value;
            //
            // Work out how many digits to converge to, normally this is
            // 2/3 of the digits in T, but if the first derivative is very
            // large convergence is slow, so we'll bump it up to full 
            // precision to prevent premature termination of the root-finding routine.
            //
            int digits = 53;
            digits /= 2;
            digits -= 1;
            if ((a < 0.125) && (Math.Abs(gamma_p_derivative(a, guess)) > 1 / Math.Sqrt(XMath.epsilon))) digits = 53;
            //
            // Go ahead and iterate:
            //
            int max_iter = max_root_iterations;
            guess = halley_iterate(new gamma_p_inverse_func(a, q, true), guess, lower, double.MaxValue, digits, max_iter);
            if (guess == lower) throw new Exception("Expected result known to be non-zero, but is smaller than the smallest available number.");
            return guess;
        }

        public static double gamma_p_inv(double a, double p)
        {
            if (a <= 0)
                throw new Exception(string.Format("Argument a in the incomplete gamma function inverse must be >= 0 (got a={0:G}).", a));
            if ((p < 0) || (p > 1))
                throw new Exception(string.Format("Probabilty must be in the range [0,1] in the incomplete gamma function inverse (got p={0:G}).", p));
            if (p == 1)
                return double.MaxValue;
            if (p == 0)
                return 0;
            bool has_10_digits;
            double guess = find_inverse_gamma(a, p, 1 - p, out has_10_digits);
            double lower = XMath.min_value;
            if (guess <= lower) guess = XMath.min_value;
            //
            // Work out how many digits to converge to, normally this is
            // 2/3 of the digits in double, but if the first derivative is very
            // large convergence is slow, so we'll bump it up to full 
            // precision to prevent premature termination of the root-finding routine.
            //
            int digits = 53;
            digits /= 2;
            digits -= 1;
            if ((a < 0.125) && (Math.Abs(gamma_p_derivative(a, guess)) > 1 / Math.Sqrt(XMath.epsilon))) digits = 51;
            //
            // Go ahead and iterate:
            //
            int max_iter = max_root_iterations;
            guess = halley_iterate(new gamma_p_inverse_func(a, p, false), guess, lower, double.MaxValue, digits, max_iter);
            if (guess == lower) throw new Exception("Expected result known to be non-zero, but is smaller than the smallest available number.");
            return guess;
        }

        public static double gamma_p_inva(double x, double p)
        {
            return gamma_inva_imp(x, p, 1 - p);
        }

        public static double gamma_q_inva(double x, double q)
        {
            return gamma_inva_imp(x, 1-q, q);
        }

        public static double gamma_p_derivative(double a, double x)
        {
            //
            // Usual error checks first:
            //
            if (a <= 0)
                throw new Exception(string.Format("Argument a to the incomplete gamma function must be greater than zero (got a={0:G}).", a));
            if (x < 0)
                throw new Exception(string.Format("Argument x to the incomplete gamma function must be >= 0 (got x={0:G}).", x));
            //
            // Now special cases:
            //
            if (x == 0)
            {
                if (a > 1) return 0.0;
                if (a == 1) return 1.0;
                throw new OverflowException();
            }
            //
            // Normal case:
            //
            double f1 = regularised_gamma_prefix(a, x);
            if (x < 1 && double.MaxValue * x < f1) throw new OverflowException();

            f1 /= x;

            return f1;
        }

        public static double lgamma(double z)
        {

            double result = 0;
            int sresult = 1;
            if (z <= 0)
            {
                // reflection formula:
                if (Math.Floor(z) == z) throw new Exception(string.Format("Evaluation of lgamma at a negative integer {0:G}.", z));

                double t = sinpx(z);
                z = -z;
                if (t < 0)
                {
                    t = -t;
                }
                else
                {
                    sresult = -sresult;
                }
                result = Math.Log(Math.PI) - lgamma(z) - Math.Log(t);
            }
            else if (z < 15)
            {
                result = lgamma_small_imp(z, z - 1, z - 2);
            }
            else if ((z >= 3) && (z < 100))
            {
                // taking the log of tgamma reduces the error, no danger of overflow here:
                result = Math.Log(gamma(z));
            }
            else
            {
                // regular evaluation:
                double zgh = z + lanczos_g - 0.5;
                result = Math.Log(zgh) - 1.0;
                result *= z - 0.5;
                result += Math.Log(lanczos_sum_expG_scaled(z));
            }
            return result;
        }

        public static double gamma_ratio(double a, double b)
        {
            return gamma_delta_ratio(a, b - a);
        }

        public static double gamma_delta_ratio(double z, double delta)
        {
            if (z <= 0)
                throw new Exception(string.Format("Gamma function ratios only implemented for positive arguments (got a={0:G}).", z));
            if (z + delta <= 0)
                throw new Exception(string.Format("Gamma function ratios only implemented for positive arguments (got b={0:G}).", z + delta));

            if (Math.Floor(delta) == delta)
            {
                if (Math.Floor(z) == z)
                {
                    //
                    // Both z and delta are integers, see if we can just use table lookup
                    // of the factorials to get the result:
                    //
                    double[] f = factorials();
                    if (z <= f.Length && (z + delta) <= f.Length)
                    {
                        return f[(int)z - 1] / f[(int)(z + delta) - 1];
                    }
                }
                if (Math.Abs(delta) < 20)
                {
                    //
                    // delta is a small integer, we can use a finite product:
                    //
                    if (delta == 0)
                        return 1;
                    if (delta < 0)
                    {
                        z -= 1;
                        double result = z;
                        while (0 != (delta += 1))
                        {
                            z -= 1;
                            result *= z;
                        }
                        return result;
                    }
                    else
                    {
                        double result = 1 / z;
                        while (0 != (delta -= 1))
                        {
                            z += 1;
                            result /= z;
                        }
                        return result;
                    }
                }
            }
            return tgamma_delta_ratio_imp_lanczos(z, delta);
        }

        public static double digamma(double x)
        {
            double result = 0;

            if (x < 0)
            {
                // Reflect:
                x = 1 - x;
                // Argument reduction for tan:
                double remainder = x - Math.Floor(x);
                if (remainder > 0.5) remainder -= 1;
                if (remainder == 0) throw new ArgumentException("Digamma function is discontinuous at negative integers");
                result = Math.PI / Math.Tan(Math.PI * remainder);
            }
            //
            // If we're above the lower-limit for the
            // asymptotic expansion then use it:
            //
            if (x >= 10) result += digamma_imp_large(x);
            else
            {
                while (x > 2)
                {
                    x -= 1;
                    result += 1 / x;
                }

                if (x < 1)
                {
                    result = -1 / x;
                    x += 1;
                }
                result += digamma_imp_1_2(x);
            }
            return result;
        }

        #region Helper functions

        private static double find_inverse_gamma(double a, double p, double q, out bool p_has_10_digits)
        {
            //
            // In order to understand what's going on here, you will
            // need to refer to:
            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //

            double result;
            p_has_10_digits = false;

            if (a == 1)
            {
                result = -Math.Log(q);
            }
            else if (a < 1)
            {
                double g = gamma(a);
                double b = q * g;
                if ((b > 0.6) || ((b >= 0.45) && (a >= 0.3)))
                {
                    // DiDonato & Morris Eq 21:
                    //
                    // There is a slight variation from DiDonato and Morris here:
                    // the first form given here is unstable when p is close to 1,
                    // making it impossible to compute the inverse of Q(a,x) for small
                    // q.  Fortunately the second form works perfectly well in this case.
                    //
                    double u;
                    if ((b * q > 1e-8) && (q > 1e-5))
                    {
                        u = Math.Pow(p * g * a, 1 / a);
                    }
                    else
                    {
                        u = Math.Exp((-q / a) - euler);
                    }
                    result = u / (1 - (u / (a + 1)));
                }
                else if ((a < 0.3) && (b >= 0.35))
                {
                    // DiDonato & Morris Eq 22:
                    double t = Math.Exp(-euler - b);
                    double u = t * Math.Exp(t);
                    result = t * Math.Exp(u);
                }
                else if ((b > 0.15) || (a >= 0.3))
                {
                    // DiDonato & Morris Eq 23:
                    double y = -Math.Log(b);
                    double u = y - (1 - a) * Math.Log(y);
                    result = y - (1 - a) * Math.Log(u) - Math.Log(1 + (1 - a) / (1 + u));
                }
                else if (b > 0.1)
                {
                    // DiDonato & Morris Eq 24:
                    double y = -Math.Log(b);
                    double u = y - (1 - a) * Math.Log(y);
                    result = y - (1 - a) * Math.Log(u) - Math.Log((u * u + 2 * (3 - a) * u + (2 - a) * (3 - a)) / (u * u + (5 - a) * u + 2));
                }
                else
                {
                    // DiDonato & Morris Eq 25:
                    double y = -Math.Log(b);
                    double c1 = (a - 1) * Math.Log(y);
                    double c1_2 = c1 * c1;
                    double c1_3 = c1_2 * c1;
                    double c1_4 = c1_2 * c1_2;
                    double a_2 = a * a;
                    double a_3 = a_2 * a;

                    double c2 = (a - 1) * (1 + c1);
                    double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                    double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
                    double c5 = (a - 1) * (-(c1_4 / 4)
                                      + (11 * a - 17) * c1_3 / 6
                                      + (-3 * a_2 + 13 * a - 13) * c1_2
                                      + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                                      + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                    double y_2 = y * y;
                    double y_3 = y_2 * y;
                    double y_4 = y_2 * y_2;
                    result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                    if (b < 1e-28) p_has_10_digits = true;
                }
            }
            else
            {
                // DiDonato and Morris Eq 31:
                double s = find_inverse_s(p, q);

                double s_2 = s * s;
                double s_3 = s_2 * s;
                double s_4 = s_2 * s_2;
                double s_5 = s_4 * s;
                double ra = Math.Sqrt(a);

                double w = a + s * ra + (s * s - 1) / 3;
                w += (s_3 - 7 * s) / (36 * ra);
                w -= (3 * s_4 + 7 * s_2 - 16) / (810 * a);
                w += (9 * s_5 + 256 * s_3 - 433 * s) / (38880 * a * ra);

                if ((a >= 500) && (Math.Abs(1 - w / a) < 1e-6))
                {
                    result = w;
                    p_has_10_digits = true;
                }
                else if (p > 0.5)
                {
                    if (w < 3 * a)
                    {
                        result = w;
                    }
                    else
                    {
                        double D = max(2.0, a * (a - 1.0));
                        double lg = lgamma(a);
                        double lb = Math.Log(q) + lg;
                        if (lb < -D * 2.3)
                        {
                            // DiDonato and Morris Eq 25:
                            double y = -lb;
                            double c1 = (a - 1) * Math.Log(y);
                            double c1_2 = c1 * c1;
                            double c1_3 = c1_2 * c1;
                            double c1_4 = c1_2 * c1_2;
                            double a_2 = a * a;
                            double a_3 = a_2 * a;

                            double c2 = (a - 1) * (1 + c1);
                            double c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
                            double c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
                            double c5 = (a - 1) * (-(c1_4 / 4)
                                              + (11 * a - 17) * c1_3 / 6
                                              + (-3 * a_2 + 13 * a - 13) * c1_2
                                              + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                                              + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

                            double y_2 = y * y;
                            double y_3 = y_2 * y;
                            double y_4 = y_2 * y_2;
                            result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
                        }
                        else
                        {
                            // DiDonato and Morris Eq 33:
                            double u = -lb + (a - 1) * Math.Log(w) - Math.Log(1 + (1 - a) / (1 + w));
                            result = -lb + (a - 1) * Math.Log(u) - Math.Log(1 + (1 - a) / (1 + u));
                        }
                    }
                }
                else
                {
                    double z = w;
                    double ap1 = a + 1;
                    if (w < 0.15 * ap1)
                    {
                        // DiDonato and Morris Eq 35:
                        double v = Math.Log(p) + lgamma(ap1);
                        double s2 = 1;
                        z = Math.Exp((v + w) / a);
                        s2 = log1p(z / ap1 * (1 + z / (a + 2)));
                        z = Math.Exp((v + z - s2) / a);
                        z = Math.Exp((v + z - s2) / a);
                        s2 = log1p(z / ap1 * (1 + z / (a + 2) * (1 + z / (a + 3))));
                        z = Math.Exp((v + z - s2) / a);
                    }

                    if ((z <= 0.01 * ap1) || (z > 0.7 * ap1))
                    {
                        result = z;
                        if (z <= 0.002 * ap1) p_has_10_digits = true;
                    }
                    else
                    {
                        // DiDonato and Morris Eq 36:
                        double ls = Math.Log(didonato_SN(a, z, 100, 1e-4));
                        double v = Math.Log(p) + lgamma(ap1);
                        z = Math.Exp((v + z - ls) / a);
                        result = z * (1 - (a * Math.Log(z) - z - v + ls) / (a - z));
                    }
                }
            }
            return result;
        }

        private static double find_inverse_s(double p, double q)
        {
            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //
            // See equation 32.
            //
            double t;
            if (p < 0.5)
            {
                t = Math.Sqrt(-2 * Math.Log(p));
            }
            else
            {
                t = Math.Sqrt(-2 * Math.Log(q));
            }
            double[] a = new double[4] { 3.31125922108741, 11.6616720288968, 4.28342155967104, 0.213623493715853 };
            double[] b = new double[5] { 1, 6.61053765625462, 6.40691597760039, 1.27364489782223, 0.3611708101884203e-1 };
            double s = t - evaluate_polynomial(a, t) / evaluate_polynomial(b, t);
            if (p < 0.5)
                s = -s;
            return s;
        }

        private static double didonato_SN(double a, double x, uint N, double tolerance)
        {
            //
            // Computation of the Incomplete Gamma Function Ratios and their Inverse
            // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
            // ACM Transactions on Mathematical Software, Vol. 12, No. 4,
            // December 1986, Pages 377-393.
            //
            // See equation 34.
            //
            double sum = 1;
            if (N >= 1)
            {
                double partial = x / (a + 1);
                sum += partial;
                for (uint i = 2; i <= N; ++i)
                {
                    partial *= x / (a + i);
                    sum += partial;
                    if (partial < tolerance)
                        break;
                }
            }
            return sum;
        }

        private class small_gamma2_series : series<double>
        {
            double result, x, apn;
            int n;
            public small_gamma2_series(double a_, double x_)
            {
                x = -x_;
                result = x;
                apn = a_ + 1;
                n = 1;
            }

            public override double next()
            {
                double r = result / (apn);
                result *= x;
                result /= ++n;
                apn += 1;
                return r;
            }
        }

        private static double finite_gamma_q(double a, double x, ref double p_derivative)
        {
            //
            // Calculates normalised Q when a is an integer:
            //
            double e = Math.Exp(-x);
            double sum = e;
            if (sum != 0)
            {
                double term = sum;
                for (uint n = 1; n < a; ++n)
                {
                    term /= n;
                    term *= x;
                    sum += term;
                }
            }
            double[] f = factorials();
            p_derivative = e * Math.Pow(x, a) / f[(int)(a - 1)];
            return sum;
        }

        private static double finite_half_gamma_q(double a, double x, ref double p_derivative)
        {
            //
            // Calculates normalised Q when a is a half-integer:
            //
            double e = erfc(Math.Sqrt(x));
            if ((e != 0) && (a > 1))
            {
                double term = Math.Exp(-x) / Math.Sqrt(Math.PI * x);
                term *= x;
                double half = 0.5;
                term /= half;
                double sum = term;
                for (uint n = 2; n < a; ++n)
                {
                    term /= n - half;
                    term *= x;
                    sum += term;
                }
                e += sum;
                p_derivative = 0;
            }
            else p_derivative = Math.Sqrt(x) * Math.Exp(-x) / XMath.root_pi;
            return e;
        }

        private static double tgamma1pm1(double dz)
        {
            double result;
            if (dz < 0)
            {
                if (dz < -0.5)
                {
                    // Best method is simply to subtract 1 from tgamma:
                    result = gamma(1 + dz) - 1;
                }
                else
                {
                    // Use expm1 on lgamma:
                    result = expm1(-log1p(dz)
                       + lgamma_small_imp(dz + 2, dz + 1, dz));
                }
            }
            else
            {
                if (dz < 2)
                {
                    // Use expm1 on lgamma:
                    result = expm1(lgamma_small_imp(dz + 1, dz, dz - 1));
                }
                else
                {
                    // Best method is simply to subtract 1 from tgamma:
                    result = gamma(1 + dz) - 1;
                }
            }

            return result;
        }

        private static double tgamma_small_upper_part(double a, double x, out double pgam, bool invert, ref double p_derivative)
        {
            //
            // Compute the full upper fraction (Q) when a is very small:
            //
            double result;
            result = tgamma1pm1(a);
            pgam = (result + 1) / a;
            double p = powm1(x, a);
            result -= p;
            result /= a;
            small_gamma2_series s = new small_gamma2_series(a, x);
            int max_iter = max_series_iterations - 10;
            p += 1;
            p_derivative = p / (pgam * Math.Exp(x));
            double init_value = invert ? pgam : 0;
            result = -p * sum_series(s, XMath.epsilon, ref max_iter, (init_value - result) / p);
            if (invert) result = -result;
            return result;
        }

        private static double igamma_temme_large(double a, double x)
        {
            double sigma = (x - a) / a;
            double phi = -log1pmx(sigma);
            double y = a * phi;
            double z = Math.Sqrt(2 * phi);
            if (x < a) z = -z;

            double[] workspace = new double[10];

            double[] C0 = 
            {
              -0.33333333333333333,
              0.083333333333333333,
              -0.014814814814814815,
              0.0011574074074074074,
              0.0003527336860670194,
              -0.00017875514403292181,
              0.39192631785224378e-4,
              -0.21854485106799922e-5,
              -0.185406221071516e-5,
              0.8296711340953086e-6,
              -0.17665952736826079e-6,
              0.67078535434014986e-8,
              0.10261809784240308e-7,
              -0.43820360184533532e-8,
              0.91476995822367902e-9,
            };
            workspace[0] = evaluate_polynomial(C0, z);

            double[] C1 = 
            {
              -0.0018518518518518519,
              -0.0034722222222222222,
              0.0026455026455026455,
              -0.00099022633744855967,
              0.00020576131687242798,
              -0.40187757201646091e-6,
              -0.18098550334489978e-4,
              0.76491609160811101e-5,
              -0.16120900894563446e-5,
              0.46471278028074343e-8,
              0.1378633446915721e-6,
              -0.5752545603517705e-7,
              0.11951628599778147e-7,
            };
            workspace[1] = evaluate_polynomial(C1, z);

            double[] C2 = 
            {
              0.0041335978835978836,
              -0.0026813271604938272,
              0.00077160493827160494,
              0.20093878600823045e-5,
              -0.00010736653226365161,
              0.52923448829120125e-4,
              -0.12760635188618728e-4,
              0.34235787340961381e-7,
              0.13721957309062933e-5,
              -0.6298992138380055e-6,
              0.14280614206064242e-6,
            };
            workspace[2] = evaluate_polynomial(C2, z);

            double[] C3 = 
            {
              0.00064943415637860082,
              0.00022947209362139918,
              -0.00046918949439525571,
              0.00026772063206283885,
              -0.75618016718839764e-4,
              -0.23965051138672967e-6,
              0.11082654115347302e-4,
              -0.56749528269915966e-5,
              0.14230900732435884e-5,
            };
            workspace[3] = evaluate_polynomial(C3, z);

            double[] C4 = 
            {
              -0.0008618882909167117,
              0.00078403922172006663,
              -0.00029907248030319018,
              -0.14638452578843418e-5,
              0.66414982154651222e-4,
              -0.39683650471794347e-4,
              0.11375726970678419e-4,
            };
            workspace[4] = evaluate_polynomial(C4, z);

            double[] C5 = 
            {
              -0.00033679855336635815,
              -0.69728137583658578e-4,
              0.00027727532449593921,
              -0.00019932570516188848,
              0.67977804779372078e-4,
              0.1419062920643967e-6,
              -0.13594048189768693e-4,
              0.80184702563342015e-5,
              -0.22914811765080952e-5,
            };
            workspace[5] = evaluate_polynomial(C5, z);

            double[] C6 = 
            {
              0.00053130793646399222,
              -0.00059216643735369388,
              0.00027087820967180448,
              0.79023532326603279e-6,
              -0.81539693675619688e-4,
              0.56116827531062497e-4,
              -0.18329116582843376e-4,
            };
            workspace[6] = evaluate_polynomial(C6, z);

            double[] C7 = 
            {
              0.00034436760689237767,
              0.51717909082605922e-4,
              -0.00033493161081142236,
              0.0002812695154763237,
              -0.00010976582244684731,
            };
            workspace[7] = evaluate_polynomial(C7, z);

            double[] C8 = 
            {
              -0.00065262391859530942,
              0.00083949872067208728,
              -0.00043829709854172101,
            };
            workspace[8] = evaluate_polynomial(C8, z);
            workspace[9] = -0.00059676129019274625;

            double result = evaluate_polynomial(workspace, 1 / a);
            result *= Math.Exp(-y) / Math.Sqrt(2 * Math.PI * a);
            if (x < a) result = -result;

            result += erfc(Math.Sqrt(y)) / 2;

            return result;
        }

        private static double gamma_incomplete_imp(double a, double x, bool normalised, bool invert, ref double p_derivative)
        {
            if (a <= 0)
                throw new Exception(string.Format("Argument a to the incomplete gamma function must be greater than zero (got a={0:G}).", a));
            if (x < 0)
                throw new Exception(string.Format("Argument x to the incomplete gamma function must be >= 0 (got x={0:G}).", x));

            double result;

            bool is_int, is_half_int;
            bool is_small_a = (a < 30) && (a <= x + 1);
            if (is_small_a)
            {
                double fa = Math.Floor(a);
                is_int = (fa == a);
                is_half_int = is_int ? false : (Math.Abs(fa - a) == 0.5);
            }
            else
            {
                is_int = is_half_int = false;
            }

            int eval_method;

            if (is_int && (x > 0.6))
            {
                // calculate Q via finite sum:
                invert = !invert;
                eval_method = 0;
            }
            else if (is_half_int && (x > 0.2))
            {
                // calculate Q via finite sum for half integer a:
                invert = !invert;
                eval_method = 1;
            }
            else if (x < 0.5)
            {
                //
                // Changeover criterion chosen to give a changeover at Q ~ 0.33
                //
                if (-0.4 / Math.Log(x) < a)
                {
                    eval_method = 2;
                }
                else
                {
                    eval_method = 3;
                }
            }
            else if (x < 1.1)
            {
                //
                // Changover here occurs when P ~ 0.75 or Q ~ 0.25:
                //
                if (x * 0.75 < a)
                {
                    eval_method = 2;
                }
                else
                {
                    eval_method = 3;
                }
            }
            else
            {
                //
                // Begin by testing whether we're in the "bad" zone
                // where the result will be near 0.5 and the usual
                // series and continued fractions are slow to converge:
                //
                bool use_temme = false;
                if (normalised && a > 20)
                {
                    double sigma = Math.Abs((x - a) / a);
                    if (a > 200)
                    {
                        //
                        // This limit is chosen so that we use Temme's expansion
                        // only if the result would be larger than about 10^-6.
                        // Below that the regular series and continued fractions
                        // converge OK, and if we use Temme's method we get increasing
                        // errors from the dominant erfc term as it's (inexact) argument
                        // increases in magnitude.
                        //
                        if (20 / a > sigma * sigma)
                            use_temme = true;
                    }
                    else
                    {
                        if (sigma < 0.4)
                            use_temme = true;
                    }
                }
                if (use_temme)
                {
                    eval_method = 5;
                }
                else
                {
                    //
                    // Regular case where the result will not be too close to 0.5.
                    //
                    // Changeover here occurs at P ~ Q ~ 0.5
                    // Note that series computation of P is about x2 faster than continued fraction
                    // calculation of Q, so try and use the CF only when really necessary, especially
                    // for small x.
                    //
                    if (x - (1 / (3 * x)) < a)
                    {
                        eval_method = 2;
                    }
                    else
                    {
                        eval_method = 4;
                        invert = !invert;
                    }
                }
            }

            switch (eval_method)
            {
                case 0:
                    {
                        result = finite_gamma_q(a, x, ref p_derivative);
                        if (normalised == false)
                            result *= gamma(a);
                        break;
                    }
                case 1:
                    {
                        result = finite_half_gamma_q(a, x, ref p_derivative);
                        if (normalised == false)
                            result *= gamma(a);
                        if (p_derivative == 0) p_derivative = regularised_gamma_prefix(a, x);
                        break;
                    }
                case 2:
                    {
                        // Compute P:
                        result = normalised ? regularised_gamma_prefix(a, x) : full_igamma_prefix(a, x);
                        p_derivative = result;
                        if (result != 0)
                        {
                            double init_value = 0;
                            if (invert)
                            {
                                init_value = -a * (normalised ? 1 : gamma(a)) / result;
                            }
                            result *= lower_gamma_series(a, x, init_value) / a;
                            if (invert)
                            {
                                invert = false;
                                result = -result;
                            }
                        }
                        break;
                    }
                case 3:
                    {
                        // Compute Q:
                        invert = !invert;
                        double g;
                        result = tgamma_small_upper_part(a, x, out g, invert, ref p_derivative);
                        invert = false;
                        if (normalised)
                            result /= g;
                        break;
                    }
                case 4:
                    {
                        // Compute Q:
                        result = normalised ? regularised_gamma_prefix(a, x) : full_igamma_prefix(a, x);
                        p_derivative = result;
                        if (result != 0)
                            result *= upper_gamma_fraction(a, x, XMath.epsilon);
                        break;
                    }
                case 5:
                    {
                        result = igamma_temme_large(a, x);
                        if (x >= a)
                            invert = !invert;
                        p_derivative = regularised_gamma_prefix(a, x);
                        break;
                    }
                default:            //needed to keep compiler happy
                    result = 0.0;
                    break;
            }

            if (normalised && (result > 1))
                result = 1;
            if (invert)
            {
                double gam = normalised ? 1 : gamma(a);
                result = gam - result;
            }
            //
            // Need to convert prefix term to derivative:
            //
            if ((x < 1) && (double.MaxValue * x < p_derivative))
            {
                // overflow, just return an arbitrarily large value:
                p_derivative = double.MaxValue / 2;
            }

            p_derivative /= x;

            return result;
        }

        class upper_incomplete_gamma_fract : series<pair<double>>
        {

            double z, a;
            int k;

            public upper_incomplete_gamma_fract(double a1, double z1)
            {
                z = z1 - a1 + 1;
                a = a1;
                k = 0;
            }

            public override pair<double> next()
            {
                ++k;
                z += 2;
                pair<double> p = new pair<double>();
                p.v1 = k * (a - k);
                p.v2 = z;
                return p;
            }
        }

        class lower_incomplete_gamma_series : series<double>
        {
            double a, z, result;
            public lower_incomplete_gamma_series(double a1, double z1)
            {
                a = a1;
                z = z1;
                result = 1;
            }

            public override double next()
            {
                double r = result;
                a += 1;
                result *= z / a;
                return r;
            }
        }

        class gamma_p_inverse_func : iterand<Tuple<double, double, double>, double>
        {
            double a, p;
            bool invert;

            public gamma_p_inverse_func(double a_, double p_, bool inv)
            {
                a = a_;
                p = p_;
                invert = inv;
                //
                // If p is too near 1 then P(x) - p suffers from cancellation
                // errors causing our root-finding algorithms to "thrash", better
                // to invert in this case and calculate Q(x) - (1-p) instead.
                //
                // Of course if p is *very* close to 1, then the answer we get will
                // be inaccurate anyway (because there's not enough information in p)
                // but at least we will converge on the (inaccurate) answer quickly.
                //
                if (p > 0.9)
                {
                    p = 1 - p;
                    invert = !invert;
                }
            }

            public override Tuple<double, double, double> next(double x)
            {
                //
                // Calculate P(x) - p and the first two derivates, or if the invert
                // flag is set, then Q(x) - q and its derivatives.
                //
                double f, f1;
                double ft = 0;
                f = gamma_incomplete_imp(a, x, true, invert, ref ft);
                f1 = ft;
                double f2;
                double div = (a - x - 1) / x;
                f2 = f1;
                if ((Math.Abs(div) > 1) && (double.MaxValue / Math.Abs(div) < f2))
                {
                    // overflow:
                    f2 = -double.MaxValue / 2;
                }
                else
                {
                    f2 *= div;
                }

                if (invert)
                {
                    f1 = -f1;
                    f2 = -f2;
                }

                return new Tuple<double, double, double>(f - p, f1, f2);
            }
        }

        private static double upper_gamma_fraction(double a, double z, double init)
        {
            upper_incomplete_gamma_fract f = new upper_incomplete_gamma_fract(a, z);
            int terms = max_series_iterations;
            return 1.0 / (z - a + 1.0 + continued_fraction_a(f, init, ref terms));
        }

        private static double lower_gamma_series(double a, double z, double init)
        {
            lower_incomplete_gamma_series s = new lower_incomplete_gamma_series(a, z);
            int terms = max_series_iterations;
            double result = sum_series(s, XMath.epsilon, ref terms, init);
            return result;
        }

        private static double lgamma_small_imp(double z, double zm1, double zm2)
        {
            double result = 0;
            if (z < XMath.epsilon)
            {
                result = -Math.Log(z);
            }
            else if ((zm1 == 0) || (zm2 == 0))
            {
                // nothing to do, result is zero....
            }
            else if (z > 2)
            {
                //
                // Begin by performing argument reduction until
                // z is in [2,3):
                //
                if (z >= 3)
                {
                    do
                    {
                        z -= 1;
                        zm2 -= 1;
                        result += Math.Log(z);
                    } while (z >= 3);
                    // Update zm2, we need it below:
                    zm2 = z - 2;
                }

                //
                // Use the following form:
                //
                // lgamma(z) = (z-2)(z+1)(Y + R(z-2))
                //
                // where R(z-2) is a rational approximation optimised for
                // low absolute error - as long as it's absolute error
                // is small compared to the constant Y - then any rounding
                // error in it's computation will get wiped out.
                //
                // R(z-2) has the following properties:
                //
                // At double: Max error found:                    4.231e-18
                // At long double: Max error found:               1.987e-21
                // Maximum Deviation Found (approximation error): 5.900e-24
                //
                double[] P = 
      {
         -0.180355685678449379109e-1,
         0.25126649619989678683e-1,
         0.494103151567532234274e-1,
         0.172491608709613993966e-1,
         -0.259453563205438108893e-3,
         -0.541009869215204396339e-3,
         -0.324588649825948492091e-4
      };
                double[] Q = 
      {
         0.1e1,
         0.196202987197795200688e1,
         0.148019669424231326694e1,
         0.541391432071720958364e0,
         0.988504251128010129477e-1,
         0.82130967464889339326e-2,
         0.224936291922115757597e-3,
         -0.223352763208617092964e-6
      };

                double Y = 0.158963680267333984375e0;

                double r = zm2 * (z + 1);
                double R = evaluate_polynomial(P, zm2);
                R /= evaluate_polynomial(Q, zm2);

                result += r * Y + r * R;
            }
            else
            {
                //
                // If z is less than 1 use recurrance to shift to
                // z in the interval [1,2]:
                //
                if (z < 1)
                {
                    result += -Math.Log(z);
                    zm2 = zm1;
                    zm1 = z;
                    z += 1;
                }
                //
                // Two approximations, on for z in [1,1.5] and
                // one for z in [1.5,2]:
                //
                if (z <= 1.5)
                {
                    //
                    // Use the following form:
                    //
                    // lgamma(z) = (z-1)(z-2)(Y + R(z-1))
                    //
                    // where R(z-1) is a rational approximation optimised for
                    // low absolute error - as long as it's absolute error
                    // is small compared to the constant Y - then any rounding
                    // error in it's computation will get wiped out.
                    //
                    // R(z-1) has the following properties:
                    //
                    // At double precision: Max error found:                1.230011e-17
                    // At 80-bit long double precision:   Max error found:  5.631355e-21
                    // Maximum Deviation Found:                             3.139e-021
                    // Expected Error Term:                                 3.139e-021

                    //
                    double Y = 0.52815341949462890625;

                    double[] P = 
         {
            0.490622454069039543534e-1,
            -0.969117530159521214579e-1,
            -0.414983358359495381969e0,
            -0.406567124211938417342e0,
            -0.158413586390692192217e0,
            -0.240149820648571559892e-1,
            -0.100346687696279557415e-2
         };
                    double[] Q = 
         {
            0.1e1,
            0.302349829846463038743e1,
            0.348739585360723852576e1,
            0.191415588274426679201e1,
            0.507137738614363510846e0,
            0.577039722690451849648e-1,
            0.195768102601107189171e-2
         };

                    double r = evaluate_polynomial(P, zm1) / evaluate_polynomial(Q, zm1);
                    double prefix = zm1 * zm2;

                    result += prefix * Y + prefix * r;
                }
                else
                {
                    //
                    // Use the following form:
                    //
                    // lgamma(z) = (2-z)(1-z)(Y + R(2-z))
                    //
                    // where R(2-z) is a rational approximation optimised for
                    // low absolute error - as long as it's absolute error
                    // is small compared to the constant Y - then any rounding
                    // error in it's computation will get wiped out.
                    //
                    // R(2-z) has the following properties:
                    //
                    // At double precision, max error found:              1.797565e-17
                    // At 80-bit long double precision, max error found:  9.306419e-21
                    // Maximum Deviation Found:                           2.151e-021
                    // Expected Error Term:                               2.150e-021
                    //
                    double Y = 0.452017307281494140625;

                    double[] P = 
         {
            -0.292329721830270012337e-1, 
            0.144216267757192309184e0,
            -0.142440390738631274135e0,
            0.542809694055053558157e-1,
            -0.850535976868336437746e-2,
            0.431171342679297331241e-3
         };
                    double[] Q = 
         {
            0.1e1,
            -0.150169356054485044494e1,
            0.846973248876495016101e0,
            -0.220095151814995745555e0,
            0.25582797155975869989e-1,
            -0.100666795539143372762e-2,
            -0.827193521891290553639e-6
         };
                    double r = zm2 * zm1;
                    double R = evaluate_polynomial(P, -zm2) / evaluate_polynomial(Q, -zm2);

                    result += r * Y + r * R;
                }
            }
            return result;
        }

        private static double tgamma_delta_ratio_imp_lanczos(double z, double delta)
        {
            double zgh = z + lanczos_g - 0.5;
            double result;
            if (Math.Abs(delta) < 10)
            {
                result = Math.Exp((0.5 - z) * log1p(delta / zgh));
            }
            else
            {
                result = Math.Pow(zgh / (zgh + delta), z - 0.5);
            }
            result *= Math.Pow(Math.E / (zgh + delta), delta);
            result *= lanczos_sum(z) / lanczos_sum(z + delta);
            return result;
        }

        private static  double full_igamma_prefix(double a, double z)
        {

            double prefix;
            double alz = a * Math.Log(z);

            if (z >= 1)
            {
                if ((alz < log_max_value) && (-z > log_min_value))
                {
                    prefix = Math.Pow(z, a) * Math.Exp(-z);
                }
                else if (a >= 1)
                {
                    prefix = Math.Pow(z / Math.Exp(z / a), a);
                }
                else
                {
                    prefix = Math.Exp(alz - z);
                }
            }
            else
            {
                if (alz > log_min_value)
                {
                    prefix = Math.Pow(z, a) * Math.Exp(-z);
                }
                else if (z / a < log_max_value)
                {
                    prefix = Math.Pow(z / Math.Exp(z / a), a);
                }
                else
                {
                    prefix = Math.Exp(alz - z);
                }
            }
            //
            // This error handling isn't very good: it happens after the fact
            // rather than before it...
            //
            if (double.IsInfinity(prefix)) throw new OverflowException("Result of incomplete gamma function is too large to represent.");

            return prefix;
        }

        private static double regularised_gamma_prefix(double a, double z)
        {
            double agh = a + lanczos_g - 0.5;
            double prefix;
            double d = ((z - a) - lanczos_g + 0.5) / agh;

            if (a < 1)
            {
                //
                // We have to treat a < 1 as a special case because our Lanczos
                // approximations are optimised against the factorials with a > 1,
                // and for high precision types especially (128-bit reals for example)
                // very small values of a can give rather eroneous results for gamma
                // unless we do this:
                //
                // TODO: is this still required?  Lanczos approx should be better now?
                //
                if (z <= log_min_value)
                {
                    // Oh dear, have to use logs, should be free of cancellation errors though:
                    return Math.Exp(a * Math.Log(z) - z - lgamma(a));
                }
                else
                {
                    // direct calculation, no danger of overflow as gamma(a) < 1/a
                    // for small a.
                    return Math.Pow(z, a) * Math.Exp(-z) / gamma(a);
                }
            }
            else if ((Math.Abs(d * d * a) <= 100) && (a > 150))
            {
                // special case for large a and a ~ z.
                prefix = a * log1pmx(d) + z * (0.5 - lanczos_g) / agh;
                prefix = Math.Exp(prefix);
            }
            else
            {
                //
                // general case.
                // direct computation is most accurate, but use various fallbacks
                // for different parts of the problem domain:
                //
                double alz = a * Math.Log(z / agh);
                double amz = a - z;
                if ((min(alz, amz) <= log_min_value) || (max(alz, amz) >= log_max_value))
                {
                    double amza = amz / a;
                    if ((min(alz, amz) / 2 > log_min_value) && (max(alz, amz) / 2 < log_max_value))
                    {
                        // compute square root of the result and then square it:
                        double sq = Math.Pow(z / agh, a / 2) * Math.Exp(amz / 2);
                        prefix = sq * sq;
                    }
                    else if ((min(alz, amz) / 4 > log_min_value) && (max(alz, amz) / 4 < log_max_value) && (z > a))
                    {
                        // compute the 4th root of the result then square it twice:
                        double sq = Math.Pow(z / agh, a / 4) * Math.Exp(amz / 4);
                        prefix = sq * sq;
                        prefix *= prefix;
                    }
                    else if ((amza > log_min_value) && (amza < log_max_value))
                    {
                        prefix = Math.Pow((z * Math.Exp(amza)) / agh, a);
                    }
                    else
                    {
                        prefix = Math.Exp(alz + amz);
                    }
                }
                else
                {
                    prefix = Math.Pow(z / agh, a) * Math.Exp(amz);
                }
            }
            prefix *= Math.Sqrt(agh / Math.E) / lanczos_sum_expG_scaled(a);
            return prefix;
        }

        class gamma_inva_t : iterand<double, double>
        {
            double z, p;
            bool invert;

            public gamma_inva_t(double z_, double p_, bool invert_)
            {
                z = z_;
                p = p_;
                invert = invert_;
            }

            public override double next(double a)
            {
                return invert ? p - gamma_q(a, z) : gamma_p(a, z) - p;
            }
        }

        private static double gamma_inva_imp(double z, double p, double q)
        {
            if (p == 0) return double.MaxValue;
            if (q == 0) return double.Epsilon;

            gamma_inva_t f = new gamma_inva_t(z, (p < q) ? p : q, (p < q) ? false : true);
            eps_tolerance tol = new eps_tolerance(53);

            double guess;
            double factor = 8;
            if (z >= 1)
            {
                guess = 1 + inverse_poisson_cornish_fisher(z, q, p);
                if (z > 5)
                {
                    if (z > 1000)
                        factor = 1.01;
                    else if (z > 50)
                        factor = 1.1;
                    else if (guess > 10)
                        factor = 1.25;
                    else
                        factor = 2;
                    if (guess < 1.1)
                        factor = 8;
                }
            }
            else if (z > 0.5) guess = z * 1.2;
            else guess = -0.4 / Math.Log(z);
            int max_iter = XMath.max_root_iterations;

            pair<double> result = bracket_and_solve_root(f, guess, factor, false, tol, ref max_iter);
            if (max_iter >= XMath.max_root_iterations)
                throw new Exception("Unable to locate the root within a reasonable number of iterations, closest approximation so far was " + result.v1.ToString());
            return (result.v1 + result.v2) / 2;
        }

        private static double digamma_imp_large(double x)
        {
            double[] P = 
            {
                0.083333333333333333333333333333333333333333333333333,
                -0.0083333333333333333333333333333333333333333333333333,
                0.003968253968253968253968253968253968253968253968254,
                -0.0041666666666666666666666666666666666666666666666667,
                0.0075757575757575757575757575757575757575757575757576,
                -0.021092796092796092796092796092796092796092796092796,
                0.083333333333333333333333333333333333333333333333333,
                -0.44325980392156862745098039215686274509803921568627
            };
            x -= 1;
            double result = Math.Log(x);
            result += 1 / (2 * x);
            double z = 1 / (x * x);
            result -= z * evaluate_polynomial(P, z);
            return result;
        }

        private static double digamma_imp_1_2(double x)
        {
            double Y = 0.99558162689208984;

            double root1 = 1569415565.0 / 1073741824;
            double root2 = (381566830.0 / 1073741824) / 1073741824;
            double root3 = 0.9016312093258695918615325266959189453125e-19;

            double[] P = 
            {    
                0.25479851061131551,
                -0.32555031186804491,
                -0.65031853770896507,
                -0.28919126444774784,
                -0.045251321448739056,
                -0.0020713321167745952
            };

            double[] Q = 
            {    
                1,
                2.0767117023730469,
                1.4606242909763515,
                0.43593529692665969,
                0.054151797245674225,
                0.0021284987017821144,
                -0.55789841321675513e-6
            };

            double g = x - root1;
            g -= root2;
            g -= root3;
            double r = evaluate_polynomial(P, x - 1) / evaluate_polynomial(Q, x - 1);
            double result = g * Y + g * r;

            return result;
        }

        #endregion
    }
}
