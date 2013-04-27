using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class weibull_distribution : distribution
    {
        double m_shape, m_scale;

        public weibull_distribution(double shape, double scale)
        {
            m_scale = scale;
            m_shape = shape;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_shape <= 0 || double.IsInfinity(m_shape)) throw new ArgumentException(string.Format("Shape argument must be a finite number > 0 (got {0:G}).", m_shape));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public override bool discrete() { return false; }

        public double shape() { return m_shape; }

        public double scale() { return m_scale; }

        public override bool LHS() { return m_shape > 1; }

        public override bool DCL() { return m_shape < 1; }

        public override bool tail_left() { return false; }

        public override bool unimodal() { return m_shape > 1; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(m_shape >= 1 ? 0 : XMath.min_value, double.MaxValue);
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double max_pdf()
        {
            if (m_shape < 1) return double.MaxValue;
            if (m_shape == 1) return 1.0 / m_scale;
            return pdf(mode());
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0)
            {
                if (m_shape < 1) return double.MaxValue;
                if (m_shape == 1) return 1;
                return 0;
            }
            if (m_shape < 1 && x <= double.Epsilon * m_scale) return double.MaxValue; 
            double result = Math.Exp(-Math.Pow(x / m_scale, m_shape));
            if (result == 0) return result;
            result *= Math.Pow(x / m_scale, m_shape) * m_shape / x;
            return result;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            double result;
            if (m_shape <= 1 && !RHS) throw new Exception("When shape parameter <= 1 there is no left-hand side for this distribution.");
            if (p == 0) result = RHS ? double.MaxValue : 0;
            else  if (m_shape == 1) result = -Math.Log(p * m_scale) * m_scale; //Exponential distribution
            else if(m_shape < 1) result = find_pdf_inv(p, double.Epsilon, quantilec(double.Epsilon), false);
            else if (RHS) result = find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            else result = find_pdf_inv(p, 0, mode(), true);
            if (double.IsInfinity(result)) result = double.MaxValue;
            return result;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return -XMath.expm1(-Math.Pow(x / m_scale, m_shape));
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return Math.Exp(-Math.Pow(x / m_scale, m_shape));
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 1) throw new OverflowException();
            return m_scale * Math.Pow(-XMath.log1p(-p), 1 / m_shape);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 1) throw new OverflowException();
            return m_scale * Math.Pow(-Math.Log(q), 1 / m_shape);
        }

        public override double mean()
        {
            return m_scale * XMath.gamma(1 + 1 / m_shape);
        }

        public override double variance()
        {
            double result = XMath.gamma(1 + 1 / m_shape);
            result += -result;
            result += XMath.gamma(1 + 2 / m_shape);
            result *= m_scale * m_scale;
            return result;
        }

        public override double mode()
        {
            if (m_shape <= 1) return 0;
            return m_scale * Math.Pow((m_shape - 1) / m_shape, 1 / m_shape);
        }

        public override double median()
        {
            return m_scale * Math.Pow(Math.Log(2.0), 1 / m_shape);
        }

        public override double skewness()
        {
            double g1 = XMath.gamma(1 + 1 / m_shape);
            double g2 = XMath.gamma(1 + 2 / m_shape);
            double g3 = XMath.gamma(1 + 3 / m_shape);
            double d = Math.Pow(g2 - g1 * g1, 1.5);
            return (2 * g1 * g1 * g1 - 3 * g1 * g2 + g3) / d;
        }

        public override double kurtosis_excess()
        {
            double g1 = XMath.gamma(1 + 1 / m_shape);
            double g2 = XMath.gamma(1 + 2 / m_shape);
            double g3 = XMath.gamma(1 + 3 / m_shape);
            double g4 = XMath.gamma(1 + 4 / m_shape);
            double g1_2 = g1 * g1;
            double g1_4 = g1_2 * g1_2;
            double d = g2 - g1_2;
            d *= d;

            double result = -6 * g1_4 + 12 * g1_2 * g2 - 3 * g2 * g2 - 4 * g1 * g3 + g4;
            result /= d;
            return result;
        }
    }
}
