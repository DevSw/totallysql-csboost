using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class binomial_distribution : distribution
    {
        ulong m_n;
        double m_p;

        public binomial_distribution(ulong n, double p)
        {
            m_n = n;
            m_p = p;
        }

        public override void check_parameters()
        {
            if (m_n <= 0 || double.IsInfinity(m_n)) throw new ArgumentException(string.Format("Number of trials must be a finite number > 0 (got {0:G}).", m_n));
            if (m_p < 0 || m_p > 1) throw new ArgumentException(string.Format("Success fraction argument must be >= 0 and <= 1 (got {0:G}).", m_p));
        }

        public override bool discrete() { return true; }

        public double success_fraction() { return m_p; }

        public ulong trials() { return m_n; }

        public enum interval_type { clopper_pearson_exact_interval = 0, jeffreys_prior_interval = 1 };

        public static double find_lower_bound_on_p(ulong trials, ulong successes, double probability, interval_type t)
        {
            if (successes == 0) return 0;

            // NOTE!!! The Clopper Pearson formula uses "successes" not
            // "successes+1" as usual to get the lower bound,
            // see http://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm
            return (t == interval_type.clopper_pearson_exact_interval) ? XMath.ibeta_inv(successes, trials - successes + 1, probability)
               : XMath.ibeta_inv(successes + 0.5, trials - successes + 0.5, probability);
        }

        public static double find_upper_bound_on_p(ulong trials, ulong successes, double probability, interval_type t)
        {
            if (trials == successes) return 1;

            // NOTE!!! The Clopper Pearson formula uses "successes" not
            // "successes+1" as usual to get the lower bound,
            // see http://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm
            return (t == interval_type.clopper_pearson_exact_interval) ? XMath.ibetac_inv(successes + 1, trials - successes, probability)
               : XMath.ibetac_inv(successes + 0.5, trials - successes + 0.5, probability);
        }

        public static ulong find_minimum_number_of_trials(ulong k, double p, double alpha)
        {
            ulong result = (ulong)XMath.ibetac_invb(k + 1, p, alpha);  // returns n - k
            return result + k;
        }

        public static ulong find_maximum_number_of_trials(ulong k, double p, double alpha)
        {
            ulong result = (ulong)XMath.ibeta_invb(k + 1, p, alpha);  // returns n - k
            return result + k;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, trials());
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, trials());
        }

//        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override double pdf(double k)
        {

            base.pdf(k);
            double n = trials();

            // Special cases of success_fraction, regardless of k successes and regardless of n trials.
            if (success_fraction() == 0)
            {  // probability of zero successes is 1:
                return k == 0.0 ? 1.0 : 0.0;
            }
            if (success_fraction() == 1)
            {  // probability of n successes is 1:
                return k == n ? 1.0 : 0.0;
            }
            // k argument may be integral, signed, or unsigned, or floating point.
            // If necessary, it has already been promoted from an integral type.
            if (n == 0) return 1; // Probability = 1 = certainty.
            if (k == 0) return Math.Pow(1 - success_fraction(), n);
            if (k == n) return Math.Pow(success_fraction(), k);  // * pow((1 - success_fraction()), (n - k)) = 1

            // Probability of getting exactly k successes
            // if C(n, k) is the binomial coefficient then:
            //
            // f(k; n,p) = C(n, k) * p^k * (1-p)^(n-k)
            //           = (n!/(k!(n-k)!)) * p^k * (1-p)^(n-k)
            //           = (tgamma(n+1) / (tgamma(k+1)*tgamma(n-k+1))) * p^k * (1-p)^(n-k)
            //           = p^k (1-p)^(n-k) / (beta(k+1, n-k+1) * (n+1))
            //           = ibeta_derivative(k+1, n-k+1, p) / (n+1)
            //
            return XMath.ibeta_derivative(k + 1, n - k + 1, success_fraction()) / (n + 1);
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            //special cases
            if (m_n == 0)
            {
                if (p == 1.0) return 0.0;
                else throw new Exception("When number of trials = 0 the inverse pdf is only defined for p = 1 (number of successes = 0)");
            }
            if (m_p == 0)
            {
                if (p == 1.0) return 0.0;
                else throw new Exception("When success fraction = 0 the inverse pdf is only defined for p = 1 (number of successes = 0)");
            }
            if (m_p == 1)
            {
                if (p == 1.0) return m_n;
                else throw new Exception("When success fraction = 1 the inverse pdf is only defined for p = 1 (number of successes = number of trials)");
            }
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            else return find_pdf_inv(p, quantile(double.Epsilon), mode(), true);
        }

        public override double cdf(double k)
        {
            base.cdf(k);
            // Cumulative Distribution Function Binomial.
            // The random variate k is the number of successes in n trials.
            // k argument may be integral, signed, or unsigned, or floating point.
            // If necessary, it has already been promoted from an integral type.

            // Returns the sum of the terms 0 through k of the Binomial Probability Density/Mass:
            //
            //   i=k
            //   --  ( n )   i      n-i
            //   >   |   |  p  (1-p)
            //   --  ( i )
            //   i=0

            // The terms are not summed directly instead
            // the incomplete beta integral is employed,
            // according to the formula:
            // P = I[1-p]( n-k, k+1).
            //   = 1 - I[p](k + 1, n - k)
            double n = trials();
            double p = success_fraction();

            if (k == n) return 1;
            // Special cases, regardless of k.
            if (p == 0)
            {  // This needs explanation:
                // the pdf is zero for all cases except when k == 0.
                // For zero p the probability of zero successes is one.
                // Therefore the cdf is always 1:
                // the probability of k or *fewer* successes is always 1
                // if there are never any successes!
                return 1;
            }
            if (p == 1)
            { // This is correct but needs explanation:
                // when k = 1
                // all the cdf and pdf values are zero *except* when k == n,
                // and that case has been handled above already.
                return 0;
            }
            //
            // P = I[1-p](n - k, k + 1)
            //   = 1 - I[p](k + 1, n - k)
            // Use of ibetac here prevents cancellation errors in calculating
            // 1-p if p is very small, perhaps smaller than machine epsilon.
            //
            // Note that we do not use a finite sum here, since the incomplete
            // beta uses a finite sum internally for integer arguments, so
            // we'll just let it take care of the necessary logic.
            //
            return XMath.ibetac(k + 1, n - k, p);
        } // binomial cdf

        public override double cdfc(double k)
        {
            base.cdfc(k);
            // Complemented Cumulative Distribution Function Binomial.
            // The random variate k is the number of successes in n trials.
            // k argument may be integral, signed, or unsigned, or floating point.
            // If necessary, it has already been promoted from an integral type.

            // Returns the sum of the terms k+1 through n of the Binomial Probability Density/Mass:
            //
            //   i=n
            //   --  ( n )   i      n-i
            //   >   |   |  p  (1-p)
            //   --  ( i )
            //   i=k+1

            // The terms are not summed directly instead
            // the incomplete beta integral is employed,
            // according to the formula:
            // Q = 1 -I[1-p]( n-k, k+1).
            //   = I[p](k + 1, n - k)
            double n = trials();
            double p = success_fraction();

            if (k == n)
            { // Probability of greater than n successes is necessarily zero:
                return 0;
            }

            // Special cases, regardless of k.
            if (p == 0)
            {
                // This need explanation: the pdf is zero for all
                // cases except when k == 0.  For zero p the probability
                // of zero successes is one.  Therefore the cdf is always
                // 1: the probability of *more than* k successes is always 0
                // if there are never any successes!
                return 0;
            }
            if (p == 1)
            {
                // This needs explanation, when p = 1
                // we always have n successes, so the probability
                // of more than k successes is 1 as long as k < n.
                // The k == n case has already been handled above.
                return 1;
            }
            //
            // Calculate cdf binomial using the incomplete beta function.
            // Q = 1 -I[1-p](n - k, k + 1)
            //   = I[p](k + 1, n - k)
            // Use of ibeta here prevents cancellation errors in calculating
            // 1-p if p is very small, perhaps smaller than machine epsilon.
            //
            // Note that we do not use a finite sum here, since the incomplete
            // beta uses a finite sum internally for integer arguments, so
            // we'll just let it take care of the necessary logic.
            //
            return XMath.ibeta(k + 1, n - k, p);
        } // binomial cdfc

        public override double quantile(double p)
        {
            base.quantile(p);
            return quantile_imp(p, 1.0 - p);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            return quantile_imp(1.0 - q, q);
        }

        public override double mean()
        {
            return trials() * success_fraction();
        }

        public override double variance()
        {
            return mean() * (1 - success_fraction());
        }

        public override double mode()
        {
            return Math.Floor(m_p * (m_n + 1));
        }

        public override double median()
        {
            return Math.Floor(m_p * m_n);
        }

        public override double skewness()
        {
            return (1.0 - 2.0 * m_p) / Math.Sqrt(m_n * m_p * (1 - m_p));
        }

        public override double kurtosis()
        {
            return 3.0 - 6.0 / m_n + 1.0 / (m_n * m_p * (1.0 - m_p));
        }

        public override double kurtosis_excess()
        {
            double q = 1.0 - m_p;
            return (1.0 - 6.0 * m_p * q) / (m_n * m_p * q);
        }

        public override double random()
        {
            return (int)Math.Round(base.random());
        }

        private double quantile_imp(double p, double q)
        { // Quantile or Percent Point Binomial function.
            // Return the number of expected successes k,
            // for a given probability p.
            //
            // Special cases:
            //
            if (p == 0)
            {  // There may actually be no answer to this question,
                // since the probability of zero successes may be non-zero,
                // but zero is the best we can do:
                return 0;
            }
            if (p == 1)
            {  // Probability of n or fewer successes is always one,
                // so n is the most sensible answer here:
                return trials();
            }
            if (p <= Math.Pow(1.0 - success_fraction(), trials()))
            { // p <= pdf(dist, 0) == cdf(dist, 0)
                return 0; // So the only reasonable result is zero.
            } // And root finder would fail otherwise.

            // Solve for quantile numerically:
            //
            double guess = XMath.inverse_binomial_cornish_fisher(trials(), success_fraction(), p, q);
            double factor = 8;
            if (trials() > 100)
                factor = 1.01; // guess is pretty accurate
            else if ((trials() > 10) && (trials() - 1 > guess) && (guess > 3))
                factor = 1.15; // less accurate but OK.
            else if (trials() < 10)
            {
                // pretty inaccurate guess in this area:
                if (guess > trials() / 64)
                {
                    guess = trials() / 4;
                    factor = 2;
                }
                else
                    guess = trials() / 1024;
            }
            else
                factor = 2; // trials largish, but in far tails.

            int max_iter = XMath.max_root_iterations;
            return inverse_discrete_quantile(p, q, guess, factor, 1.0, max_iter);
        } // quantile_imp

    }
}
