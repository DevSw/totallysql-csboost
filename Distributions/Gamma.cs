using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class gamma_distribution : distribution
    {
        double m_shape;
        double m_scale;

        public gamma_distribution(double shape, double scale)
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

        public double shape()
        {
            return m_shape;
        }

        public double scale()
        {
            return m_scale;
        }

        public override bool LHS() { return m_shape > 1; }

        public override bool unimodal()
        {
            return m_shape > 1;
        }

        public override bool DCL()
        {
            return m_shape < 1;
        }

        public override bool tail_left() { return false; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(XMath.min_value, double.MaxValue);
        }

        public override double max_pdf()
        {
            if (m_shape < 1) return double.MaxValue;
            if (m_shape == 1) return 1 / m_scale;
            return pdf(mode());
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0) return 0;
            return XMath.gamma_p_derivative(m_shape, x / m_scale) / m_scale;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            if (!RHS && !LHS()) throw new Exception("The Gamma distribution has no LHS when the shape parameter is <= 1 - so no LHS inverse PDF is available");
            base.pdf_inv(p, RHS);
            if (m_shape <= 1)
            {
                if (p >= double.MaxValue) return 0;
                return find_pdf_inv(p, double.Epsilon, quantilec(double.Epsilon), false);
            }
            else if (RHS) return find_pdf_inv(p, mode(), quantilec(double.Epsilon), false);
            else return find_pdf_inv(p, 0, mode(), true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            return XMath.gamma_p(m_shape, x / m_scale);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return XMath.gamma_q(m_shape, x / m_scale);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 1) throw new OverflowException();
            return XMath.gamma_p_inv(m_shape, p) * m_scale;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            if (q == 0) throw new OverflowException();
            return XMath.gamma_q_inv(m_shape, q) * m_scale;
        }

        public override double mean()
        {
            return m_shape * m_scale;
        }

        public override double variance()
        {
            return m_shape * m_scale * m_scale;
        }

        public override double mode()
        {
            if (m_shape < 1)
                throw new Exception(string.Format("gamma.mode: The mode of the gamma distribution is only defined for values of the shape parameter >= 1, but got {0:G}.", m_shape));
            return (m_shape - 1) * m_scale;
        }

        //Median supplied by base class

        public override double skewness()
        {
            return 2.0 / Math.Sqrt(m_shape);
        }

        public override double kurtosis_excess()
        {
            return 6.0 / m_shape;
        }

        //Kurtosis supplied by base class
    }
}
