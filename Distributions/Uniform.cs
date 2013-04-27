using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class uniform_distribution : distribution
    {
        double m_lower, m_upper;

        public uniform_distribution(double lower, double upper)
        {
            m_lower = lower;
            m_upper = upper;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_lower)) throw new ArgumentException(string.Format("Uniform Distribution: Both bounds must be finite numbers (got {0:G}).", m_lower));
            if (double.IsInfinity(m_upper)) throw new ArgumentException(string.Format("Uniform Distribution: Both bounds must be finite numbers (got {0:G}).", m_upper));
            if (m_upper <= m_lower) throw new ArgumentException(string.Format("Uniform Distribution: Upper bound must be > lower bound (got lower bound = {0:G}, upper bound = {1:G}).", m_lower, m_upper));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return true; }

        public double lower() { return m_lower; }

        public double upper() { return m_upper; }

        public override bool RHS() { return false; }

        public override bool LHS() { return false; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override bool unimodal() { return false; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(m_lower, m_upper);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x < m_lower || x > m_upper) return 0;
            return 1.0 / (m_upper - m_lower);
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            throw new Exception("Uniform Distribution: inverse pdf is not defined");
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if ((x <= m_lower)) return 0;
            if (x >= m_upper) return 1;
            return (x - m_lower) / (m_upper - m_lower);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x <= m_lower) return 1;
            if (x >= m_upper) return 0;
            return (m_upper - x) / (m_upper - m_lower);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return m_lower;
            if (p == 1) return m_upper;
            return p * (m_upper - m_lower) + m_lower;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 0) return m_upper;
            if (q == 1) return m_lower;
            return q * (m_lower - m_upper) + m_upper;
        }

        public override double mean()
        {
            return (m_lower + m_upper) / 2;
        }

        public override double variance()
        {
            return (m_upper - m_lower) * (m_upper - m_lower) / 12;
        }

        public override double mode()
        {
            return m_lower;
        }

        public override double median()
        {
            return mean();
        }

        public override double skewness()
        {
            return 0;
        }

        public override double kurtosis()
        {
            return 3.0 - 6.0 / 5;
        }

        public override double kurtosis_excess()
        {
            return 6.0 / 5;
        }
    }
}
