using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class triangular_distribution : distribution
    {
        double m_lower, m_mode, m_upper;

        public triangular_distribution(double lower, double mode, double upper)
        {
            m_lower = lower;
            m_mode = mode;
            m_upper = upper;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_lower)) throw new ArgumentException(string.Format("All arguments must be a finite number (got {0:G}).", m_lower));
            if (double.IsInfinity(m_mode)) throw new ArgumentException(string.Format("All arguments must be a finite number (got {0:G}).", m_mode));
            if (double.IsInfinity(m_upper)) throw new ArgumentException(string.Format("All arguments must be a finite number (got {0:G}).", m_upper));
            if (m_upper <= m_lower) throw new ArgumentException(string.Format("Upper bound must be > lower bound of triangle (got lower bound = {0:G}, upper bound = {1:G}).", m_lower, m_upper));
            if (m_mode < m_lower) throw new ArgumentException(string.Format("Mode argument must be >= lower bound of triangle (got lower bound = {0:G}, mode = {1:G}).", m_lower, m_mode));
            if (m_mode > m_upper) throw new ArgumentException(string.Format("Mode argument must be <= upper bound of triangle (got upper bound = {0:G}, mode = {1:G}).", m_upper, m_mode));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return (m_mode - m_lower == m_upper - m_mode); }

        public double lower() { return m_lower; }

        public override double mode() { return m_mode; }

        public double upper() { return m_upper; }

        public override bool RHS() { return m_mode < m_upper; }

        public override bool LHS() { return m_mode > m_lower; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

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
            if (x == m_lower) return m_mode == m_lower ? 2 / (m_upper - m_lower) : 0.0;
            if (x == m_upper) return m_mode == m_upper ? 2 / (m_upper - m_lower) : 0.0;
            if (x <= m_mode) return 2 * (x - m_lower) / ((m_upper - m_lower) * (m_mode - m_lower));
            return 2 * (m_upper - x) / ((m_upper - m_lower) * (m_upper - m_mode));
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? m_upper : m_lower;
            if (RHS && m_mode == m_upper) throw new Exception("This distribution has no RHS when mode = upper bound.");
            if (!RHS && m_mode == m_lower) throw new Exception("This distribution has no LHS when mode = lower bound");
            if (RHS) return (m_mode - m_upper) * (p / pdf(m_mode)) + m_upper;
            return (m_mode - m_lower) * (p / pdf(m_mode)) + m_lower;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if ((x <= m_lower)) return 0;
            if (x >= m_upper) return 1;
            if (x <= m_mode) return ((x - m_lower) * (x - m_lower)) / ((m_upper - m_lower) * (m_mode - m_lower));
            return 1 - (m_upper - x) * (m_upper - x) / ((m_upper - m_lower) * (m_upper - m_mode));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x <= m_lower) return 1;
            if (x >= m_upper) return 0;
            if (x <= m_mode) return 1 - ((x - m_lower) * (x - m_lower)) / ((m_upper - m_lower) * (m_mode - m_lower));
            return (m_upper - x) * (m_upper - x) / ((m_upper - m_lower) * (m_upper - m_mode));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return m_lower;
            if (p == 1) return m_upper;
            double p0 = (m_mode - m_lower) / (m_upper - m_lower);
            double q = 1 - p;
            if (p < p0) return Math.Sqrt((m_upper - m_lower) * (m_mode - m_lower) * p) + m_lower;
            if (p == p0) return m_mode;
            return m_upper - Math.Sqrt((m_upper - m_lower) * (m_upper - m_mode) * q);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);

            double p = 1 - q;
            double p0 = (m_mode - m_lower) / (m_upper - m_lower);
            if (p < p0)
            {
                double s = (m_upper - m_lower) * (m_mode - m_lower);
                s *= p;
                return Math.Sqrt((m_upper - m_lower) * (m_mode - m_lower) * p) + m_lower;
            }
            if (p == p0) return m_mode;
            return m_upper - Math.Sqrt((m_upper - m_lower) * (m_upper - m_mode) * q);
        }

        public override double mean()
        {
            return (m_lower + m_mode + m_upper) / 3;
        }

        public override double variance()
        {
            return (m_lower * m_lower + m_upper * m_upper + m_mode * m_mode - m_lower * m_upper - m_lower * m_mode - m_upper * m_mode) / 18;
        }

        public override double median()
        {
            if (m_mode < (m_upper - m_lower) / 2) return m_lower + Math.Sqrt((m_upper - m_lower) * (m_mode - m_lower)) / XMath.root_two;
            return m_upper - Math.Sqrt((m_upper - m_lower) * (m_upper - m_mode)) / XMath.root_two;
        }

        public override double skewness()
        {
            return XMath.root_two * (m_lower + m_upper - 2 * m_mode) * (2 * m_lower - m_upper - m_mode) * (m_lower - 2 * m_upper + m_mode) /
            (5 * Math.Pow((m_lower * m_lower + m_upper + m_upper + m_mode * m_mode - m_lower * m_upper - m_lower * m_mode - m_upper * m_mode), 3.0 / 2.0));
        }

        public override double kurtosis()
        {
            return 12.0 / 5;
        }

        public override double kurtosis_excess()
        {
            return -3.0 / 5;
        }
    }
}
