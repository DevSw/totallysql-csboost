using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class beta_distribution : distribution
    {
        double m_alpha;
        double m_beta;

        public beta_distribution(double alpha, double beta)
        {
            m_alpha = alpha;
            m_beta = beta;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_alpha <= 0 || double.IsInfinity(m_alpha)) throw new ArgumentException(string.Format("Alpha argument must be a finite number > 0 (got {0:G}).", m_alpha));
            if (m_beta <= 0 || double.IsInfinity(m_beta)) throw new ArgumentException(string.Format("Beta argument must be a finite number > 0 (got {0:G}).", m_beta));
        }

        public override bool discrete() { return false; }

        public override bool inverted()
        {
            return m_alpha < 1 && m_beta < 1;
        }

        public override bool unimodal()
        {
            return (m_alpha >= 1 && m_beta >= 1);
        }

        public double alpha()
        {
            return m_alpha;
        }

        public double beta()
        {
            return m_beta;
        }

        public override bool LHS() { return !(m_alpha < 1 && m_beta >= 1 || m_alpha == 1 && m_beta > 1); }

        public override bool RHS() { return !(m_beta < 1 && m_alpha >= 1 || m_beta == 1 && m_alpha > 1); }

        public override bool DCL() { return m_alpha < 1; }

        public override bool DCR() { return m_beta < 1; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override bool symmetric() { return m_alpha == m_beta; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, 1);
        }

        public override XMath.pair<double> support()
        {
            double lo, hi;
            lo = m_alpha < 1 ? XMath.epsilon : 0;
            hi = m_beta < 1 ? 1 - XMath.epsilon : 1;
            return new XMath.pair<double>(lo, hi);
        }

        public override double pdf(double x)
        {
            base.pdf(x);
            return XMath.ibeta_derivative(m_alpha, m_beta, x);
        }

        public override double max_pdf()
        {
            if (m_alpha < 1 || m_beta < 1) return double.PositiveInfinity;
            if (m_alpha > 1 && m_beta > 1) return pdf(mode());
            if (m_alpha == 1 && m_beta == 1) return 1;
            if (m_alpha < 1 || (m_alpha == 1 && m_beta >= 1)) return pdf(0);
            return pdf(1);
        }

        public override double pdf_inv(double p, bool rhs)
        {
            if(rhs && !RHS()) throw new Exception(string.Format("The beta distribution has no RHS when beta {0:G} 1 and alpha {1:G} 1.", m_beta < 1 ? "<" : "==", m_beta < 1 ? ">=" : ">"));
            if (!rhs && !LHS()) throw new Exception(string.Format("The beta distribution has no LHS when alpha {0:G} 1 and beta {1:G} 1.", m_alpha < 1 ? "<" : "==", m_alpha< 1 ? ">=" : ">"));
            base.pdf_inv(p, rhs);
            double result = 0, midpoint;
            if(m_alpha == 1 && m_beta == 1) throw new Exception("Inverse pdf is not defined for this distribution when alpha = beta = 1");
            if (m_alpha < 1 && m_beta < 1)
            {
                midpoint = antimode();
                double pmin = pdf(midpoint);
                if (p < pmin) throw new ArgumentException(string.Format("Minimum pdf value for this distribution is {0:G} (got {1:G}).", pmin, p));
                if (rhs) result = find_pdf_inv(p, midpoint, support().v2, true);
                else result = find_pdf_inv(p, support().v1, midpoint, false);
            }
            else if (m_alpha > 1 && m_beta > 1)
            {
                midpoint = mode();
                if (rhs) result = find_pdf_inv(p, midpoint, support().v2, false);
                else result = find_pdf_inv(p, support().v1, midpoint, true);
            }
            else if (m_alpha < 1 || (m_alpha == 1 && m_beta >= 1)) result = find_pdf_inv(p, support().v1, support().v2, false);
            else result = find_pdf_inv(p, support().v1, support().v2, true);
            return result;
        }

        public override double cdf(double x)
        {
            base.cdf(x);
            if (x == 0) return 0;
            else if (x == 1) return 1;
            return XMath.ibeta(m_alpha, m_beta, x);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            if (x == 0) return 1;
            else if (x == 1) return 0;
            return XMath.ibetac(m_alpha, m_beta, x);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            if (p == 0) return 0;
            if (p == 1) return 1;
            return XMath.ibeta_inv(m_alpha, m_beta, p);
        }

        public override double quantilec(double p)
        {
            base.quantilec(p);
            if (p == 0) return 1;
            if (p == 1) return 0;
            return XMath.ibetac_inv(m_alpha, m_beta, p);
        }

        public override double mean()
        {
            return m_alpha / (m_alpha + m_beta);
        }

        public override double variance()
        {
            double a = m_alpha, b = m_beta;
            return (a * b) / ((a + b) * (a + b) * (a + b + 1));
        }

        public override double antimode()   //Same as mode but represents lowest point in U-shaped distribution.
        {
            if (m_alpha >= 1) throw new Exception(string.Format("beta.mode: antimode undefined for alpha = {0:G}, must be < 1!", m_alpha));
            if (m_beta >= 1) throw new Exception(string.Format("beta.mode: antimode undefined for beta = {0:G}, must be < 1!", m_beta));
            double a = m_alpha, b = m_beta;
            return (a - 1) / (a + b - 2);
        }

        public override double mode()
        {
            if (m_alpha < 1) throw new Exception(string.Format("beta.mode: mode undefined for alpha = {0:G}, must be >= 1!", m_alpha));
            if (m_beta < 1) throw new Exception(string.Format("beta.mode: mode undefined for beta = {0:G}, must be >= 1!", m_beta));

            double a = m_alpha, b = m_beta;
            return (a - 1) / (a + b - 2);
        } // mode

        //median supplied by base class

        public override double skewness()
        {
            double a = m_alpha, b = m_beta;
            return (2 * (b - a) * Math.Sqrt(a + b + 1)) / ((a + b + 2) * Math.Sqrt(a * b));
        }

        public override double kurtosis_excess()
        {
            double a = m_alpha, b = m_beta;
            double a_2 = a * a;
            double n = 6 * (a_2 * a - a_2 * (2 * b - 1) + b * b * (b + 1) - 2 * a * b * (b + 2));
            double d = a * b * (a + b + 2) * (a + b + 3);
            return n / d;
        }

        //kurtosis supplied by base class

        public override void setup_bars()
        {
            if (m_alpha == 1 && m_beta == 1) return;                                //Uniform distibution
            base.setup_bars();
        }

        public override double random()
        {
            if (m_alpha == 1 && m_beta == 1) return distribution.rand.NextDouble();             //Uniform distibution
            if (m_alpha > 1 || m_beta > 1) return base.random();
            return(base.random());
        }
    }
}
