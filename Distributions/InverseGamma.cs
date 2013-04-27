using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class inverse_gamma_distribution : distribution
    {
        double m_shape, m_scale;

        public inverse_gamma_distribution(double shape, double scale)
        {
            m_shape = shape;
            m_scale = scale;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_shape <= 0 || double.IsInfinity(m_shape)) throw new ArgumentException(string.Format("Inverse Gamma Distribution: Shape argument must be a finite number > 0 (got {0:G}).", m_shape));
            if (m_scale <= 0 || double.IsInfinity(m_scale)) throw new ArgumentException(string.Format("Inverse Gamma Distribution: Scale argument must be a finite number > 0 (got {0:G}).", m_scale));
        }

        public override bool discrete() { return false; }

        public double shape() { return m_shape; }

        public double scale() { return m_scale; }

        public override bool LHS() { return true; }

        public override bool tail_left() { return false; }

        public override bool unimodal()
        {
            return true;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            if (x == 0) return 0;

            double result = m_scale / x;
            if(result < double.Epsilon) return 0;  
            result = XMath.gamma_p_derivative(m_shape, result) * m_scale;
            if (result != 0)
            {
                if (x < 0)
                {
                    double lim = double.MaxValue * x;
                    if (lim < result) throw new OverflowException();
                    result /= x;
                    if (lim < result) throw new OverflowException();
                    result /= x;
                }
                result /= (x * x);
            }
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
            if (p == max_pdf()) return mode();

            double ubound = RHS ? double.MaxValue : mode();
            double lbound = RHS ? mode() : 0;
            return find_pdf_inv(p, lbound, ubound, RHS ? false : true);
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x == 0) return 0;
            return XMath.gamma_q(m_shape, m_scale / x);        
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x == 0) return 1;
            return XMath.gamma_p(m_shape, m_scale / x); 
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double result = XMath.gamma_q_inv(m_shape, p);
            if(result == 0) return double.PositiveInfinity;
            result = m_scale / result;
            return result;
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double result = XMath.gamma_p_inv(m_shape, q);
            if (result == 0) return double.PositiveInfinity;
            result = m_scale / result;
            return result;
        }

        public override double mean()
        {
            if (m_shape <= 1) throw new Exception(string.Format("Inverse Gamma Distribution: The mean is only defined when shape > 1 (got {0:G}).", m_shape));
            return m_scale / (m_shape - 1.0);
        }

        public override double variance()
        {
            if (m_shape <= 2) throw new Exception(string.Format("Inverse Gamma Distribution: The variance is only defined when shape > 2 (got {0:G}).", m_shape));
            return (m_scale * m_scale) / ((m_shape - 1) * (m_shape - 1) * (m_shape - 2));
        }

        public override double mode()
        {
            return m_scale / (m_shape + 1);
        }

        //Median supplied by base class

        public override double skewness()
        {
            if (m_shape <= 3) throw new Exception(string.Format("Inverse Gamma Distribution: The skewness is only defined when shape > 3 (got {0:G}).", m_shape));
            return (4 * Math.Sqrt(m_shape - 2)) / (m_shape - 3);
        }

        //Kurtosis supplied by base class

        public override double kurtosis_excess()
        {
            if (m_shape <= 4) throw new Exception(string.Format("Inverse Gamma Distribution: The kurtosis is only defined when shape > 4 (got {0:G}).", m_shape));
            return (30 * m_shape - 66) / ((m_shape - 3) * (m_shape - 4));
        }
    }
}
