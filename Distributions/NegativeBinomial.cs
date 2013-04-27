using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class negative_binomial_distribution : distribution
    {
        double m_r, m_p;

        public negative_binomial_distribution(double r, double p)
        {
            m_r = r;
            m_p = p;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_r <= 0 || double.IsInfinity(m_r)) throw new ArgumentException(string.Format("Number of successes must be a finite number > 0 (got {0:G}).", m_r));
            if (m_p < 0 || m_p > 1) throw new ArgumentException(string.Format("Success fraction argument must be >= 0 and <= 1 (got {0:G}).", m_p));
        }

        public double success_fraction() { return m_p; }

        public double successes() { return m_r; }

        public override bool discrete() { return true; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public double find_lower_bound_on_p(double trials, double successes, double alpha)
        {
            double failures = trials - successes;
            return XMath.ibeta_inv(successes, failures + 1, alpha);
        }

        public double find_upper_bound_on_p(double trials, double successes, double alpha)
        {
            double failures = trials - successes;
            return XMath.ibetac_inv(successes, failures, alpha);
        }

        public double find_minimum_number_of_trials(double k, double p, double alpha)
        {
            return XMath.ibeta_inva(k + 1, p, alpha) + k;
        }

        public double find_maximum_number_of_trials(double k, double p, double alpha)
        {
            return XMath.ibetac_inva(k + 1, p, alpha) + k;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double mean()
        {
            return m_r * (1 - m_p) / m_p;
        }

        public override double variance()
        {
            return m_r * (1 - m_p) / (m_p * m_p);
        }

        public override double mode()
        {
            return Math.Floor((m_r - 1) * (1 - m_p) / m_p);
        }

        //Median supplied by base class

        public override double skewness()
        {
            return (2 - m_p) / Math.Sqrt(m_r * (1 - m_p));
        }

        public override double kurtosis()
        {
            return 3 + (6 / m_r) + ((m_p * m_p) / (m_r * (1 - m_p)));
        }

        public override double kurtosis_excess()
        {
            return (6 - m_p * (6 - m_p)) / (m_r * (1 - m_p));
        }

        public override double pdf(double k)
        {
            base.pdf(k);
            return (m_p / (m_r + k)) * XMath.ibeta_derivative(m_r, k + 1, m_p);
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            return find_pdf_inv(p, 0, mode(), true);
        }

        public override double cdf(double k)
        {
            base.cdf(k);
            return XMath.ibeta(m_r, k + 1,m_p);
        }

        public override double cdfc(double k)
        {
            base.cdfc(k);
            return XMath.ibetac(m_r, k + 1, m_p);
        }

        public override double quantile(double P)
        {
            base.quantile(P);
            double p = m_p;
            double r = m_r;

            if (P == 1) throw new OverflowException();
            if (P == 0) return 0; // Total trials will be just dist.successes.
            if (P <= Math.Pow(p, r)) return 0;

            double guess = 0;
            double factor = 5;
            if (r * r * r * P * p > 0.005)
                guess = XMath.inverse_negative_binomial_cornish_fisher(r, p, 1.0 - p, P, 1.0 - P);

            if (guess < 10) guess = Math.Min(r * 2.0, 10.0);
            else factor = (1 - P < Math.Sqrt(XMath.epsilon)) ? 2.0 : (guess < 20 ? 1.2 : 1.1);
            int max_iter = XMath.max_root_iterations;
            return inverse_discrete_quantile(P, 1 - P, guess, factor, 1.0, max_iter);
        }

        public override double quantilec(double Q)
        {
            base.quantilec(Q);
            double p = m_p;
            double r = m_r;
            if (Q == 1) return 0;
            if (-Q <= XMath.powm1(p, r)) return 0;
            if (Q == 0) throw new OverflowException();
            double guess = 0;
            double factor = 5;
            if (r * r * r * (1 - Q) * p > 0.005) guess = XMath.inverse_negative_binomial_cornish_fisher(r, p, 1.0 - p, 1.0 - Q, Q);

            if (guess < 10) guess = Math.Min(r * 2, 10.0);
            else factor = (Q < Math.Sqrt(XMath.epsilon)) ? 2.0 : (guess < 20 ? 1.2 : 1.1);
            int max_iter = XMath.max_root_iterations;
            return inverse_discrete_quantile(1 - Q, Q, guess, factor, 1.0, max_iter);
        }

        public override double random()
        {
            return Math.Round(base.random());
        }
    }
}
