using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class laplace_distribution : distribution
    {
        double m_location, m_scale;

        public laplace_distribution(double location, double scale)
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

        public override bool symmetric() { return true; }

        public double location() { return m_location; }

        public double scale() { return m_scale; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return range();
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            double exponent = x - m_location;
            if (exponent > 0) exponent = -exponent;
            exponent /= m_scale;

            double result = Math.Exp(exponent);
            result /= 2 * m_scale;

            return result;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            double x = RHS ? m_location + Math.Abs(m_scale * Math.Log(p * 2 * m_scale)) : m_location - Math.Abs(m_scale * Math.Log(p * 2 * m_scale));
            return x;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if(double.IsNegativeInfinity(x)) return 0;
            if(double.IsPositiveInfinity(x)) return 1;
            double result;
            if (x < m_location)
                result = Math.Exp((x - m_location) / m_scale) / 2;
            else
                result = 1 - Math.Exp((m_location - x) / m_scale) / 2;

            return result;
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (double.IsNegativeInfinity(x)) return 1;
            if (double.IsPositiveInfinity(x)) return 0;
            double result;
            if (-x < m_location)
                result = Math.Exp((-x - m_location) / m_scale) / 2;
            else
                result = 1 - Math.Exp((m_location + x) / m_scale) / 2;

            return result;
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double result;
            if (p == 0) return double.NegativeInfinity;
            if (p == 1) return double.PositiveInfinity;
            if (p - 0.5 < 0.0)
                result = m_location + m_scale * Math.Log(p * 2);
            else
                result = m_location - m_scale * Math.Log(-p * 2 + 2);

            return result;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double result;
            if (q == 0) return double.PositiveInfinity;
            if (q == 1) return double.NegativeInfinity;
            if (0.5 - q < 0.0)
                result = m_location + m_scale * Math.Log(-q * 2 + 2);
            else
                result = m_location - m_scale * Math.Log(q * 2);

            return result;
        }

        public override double mean()
        {
            return m_location;
        }

        public override double standard_deviation()
        {
            return XMath.root_two * m_scale;
        }

        public override double mode()
        {
            return m_location;
        }

        public override double median()
        {
            return m_location;
        }

        public override double skewness()
        {
            return 0;
        }

        public override double kurtosis()
        {
            return 6;
        }

        public override double kurtosis_excess()
        {
            return 3;
        }
    }
}
