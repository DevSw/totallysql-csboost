using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class lognormal_distribution : distribution
    {
        double m_location, m_scale;

        public lognormal_distribution(double location, double scale)
        {
            m_location = location;
            m_scale = scale;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (double.IsInfinity(m_location)) throw new ArgumentException(string.Format("Location argument must be a finite number (got {0:G}).", m_location));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public override bool discrete() { return false; }

        public double location() { return m_location; }

        public double scale() { return m_scale; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override bool tail_left() { return false; }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0) return 0;

            double exponent = Math.Log(x) - m_location;
            exponent *= -exponent;
            exponent /= 2 * m_scale * m_scale;
            double result = Math.Exp(exponent);
            result /= m_scale * XMath.root_two_pi * x;
            return result;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : 0;
            if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            return find_pdf_inv(p, 0, mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x == 0) return 0;
            normal_distribution norm = new normal_distribution(m_location, m_scale);
            return norm.cdf(Math.Log(x));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x == 0) return 0;
            normal_distribution norm = new normal_distribution(m_location, m_scale);
            return norm.cdfc(Math.Log(x));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return 0;
            if (p == 1) throw new OverflowException();
            normal_distribution norm = new normal_distribution(m_location, m_scale);
            return Math.Exp(norm.quantile(p));
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 1) return 0;
            if (q == 0) throw new OverflowException();
            normal_distribution norm = new normal_distribution(m_location, m_scale);
            double result = Math.Exp(norm.quantilec(q));
            if (double.IsPositiveInfinity(result)) result = double.MaxValue;
            else if (double.IsNegativeInfinity(result)) result = -double.MaxValue;
            return result;
        }

        public override double mean()
        {
            return Math.Exp(m_location + m_scale * m_scale / 2.0);
        }

        public override double variance()
        {
            return XMath.expm1(m_scale * m_scale) * Math.Exp(2 * m_location + m_scale * m_scale);
        }

        public override double mode()
        {
            return Math.Exp(m_location - m_scale * m_scale);
        }

        public override double median()
        {
            return Math.Exp(m_location);
        }

        public override double skewness()
        {
            double ss = m_scale * m_scale;
            double ess = Math.Exp(ss);
            return (ess + 2) * Math.Sqrt(XMath.expm1(ss));
        }

        public override double kurtosis()
        {
            double ss = m_scale * m_scale;
            return Math.Exp(4 * ss) + 2 * Math.Exp(3 * ss) + 3 * Math.Exp(2 * ss) - 3;
        }

        public override double kurtosis_excess()
        {
            double ss = m_scale * m_scale;
            return Math.Exp(4 * ss) + 2 * Math.Exp(3 * ss) + 3 * Math.Exp(2 * ss) - 6;
        }
    }
}
