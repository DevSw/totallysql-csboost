using System;
using System.Collections.Generic;

using System.Text;
using System.Runtime.Serialization;
using Microsoft.SqlServer.Server;

namespace CSBoost.Distributions
{
    public abstract class distribution : IBinarySerialize
    {
        protected static readonly Random rand = new Random();
        double ymax = -1;

        protected virtual int barcount { get { return 256; } }
        public abstract void check_parameters();
        public abstract bool discrete();
        public virtual bool inverted() { return false; }
        public virtual double antimode()
        {
            throw new NotImplementedException();
        }
        public virtual bool unimodal() { return true; }
        public virtual bool RHS() { return true; }          //Does the distribution have a right-hand side
        public virtual bool LHS() { return true; }          //Does the distribution have a left-hand side
        public virtual bool DCL() { return false; }          //Does the distribution have a vertical discontinuity on the left
        public virtual bool DCR() { return false; }         //Does the distribution have a vertical discontinuity on the right
        public virtual bool strictly_increasing() { return !unimodal() && !inverted() && LHS(); }
        public virtual bool strictly_decreasing() { return !unimodal() && !inverted() && RHS(); }
        public virtual bool tail_right() { return true; }
        public virtual bool tail_left() { return true; }
        public virtual bool symmetric() { return false; }
        public abstract XMath.pair<double> range();
        public abstract XMath.pair<double> support();
        public virtual double pdf(double x)
        {
            XMath.pair<double> r = range();
            if (x < r.v1 || x > r.v2)
                throw new ArgumentException(string.Format("{0:G}.pdf: Permitted values for x are {1:G} to {2:G} (got {3:G}).", this.GetType().Name, r.v1, r.v2, x));
            return 0;
        }
        public virtual double max_pdf() 
        {
            if(unimodal()) return pdf(mode());
            if (DCL() || DCR()) return double.MaxValue;
            return Math.Max(pdf(support().v1), pdf(support().v2));
        }
        public virtual double min_pdf()
        {
            if (inverted()) return pdf(antimode());
            return Math.Min(pdf(support().v1), pdf(support().v2));
        }
        public virtual double pdf_inv(double p, bool RHS)
        {
            if (ymax == -1) ymax = max_pdf();
            if (p < 0) throw new ArgumentException(string.Format("Argument must be > 0 (got {0:G}", p));
            if (p > ymax) throw new ArgumentException(string.Format("Argument too high: the maximum pdf value for this disrtibution is {0:G}", ymax));
            return double.NaN;
        }
        public virtual double cdf(double x)
        {
            XMath.pair<double> r = range();
            if (x < r.v1 || x > r.v2)
                throw new ArgumentException(string.Format("{0:G}.cdf: Permitted values for x are {1:G} to {2:G} (got {3:G}).", this.GetType().Name, r.v1, r.v2, x));
            return 0;
        }
        public virtual double cdfc(double x)
        {
            XMath.pair<double> r = range();
            if (x < r.v1 || x > r.v2)
                throw new ArgumentException(string.Format("{0:G}.cdfc: Permitted values for x are {1:G} to {2:G} (got {3:G}).", this.GetType().Name, r.v1, r.v2, x));
            return 0;
        }
        public virtual double quantile(double p)
        {
            if (p < 0 || p > 1)
                throw new ArgumentException(string.Format("{0:G}.quantile: Permitted values for p are {1:G} to {2:G} (got {3:G}).", this.GetType().Name, 0, 1, p));
            return 0;
        }
        public virtual double quantilec(double q)
        {
            if (q < 0 || q > 1)
                throw new ArgumentException(string.Format("{0:G}.quantile: Permitted values for q are {1:G} to {2:G} (got {3:G}).", this.GetType().Name, 0, 1, q));
            return 0;
        }
        public abstract double mean();
        public virtual double variance() { double s = standard_deviation(); return s * s; }  //either this or standard_deviation must be over-ridden to avoid endless loop
        public abstract double mode();
        public virtual double median() { return quantile(0.5); }
        public abstract double skewness();
        public virtual double kurtosis() { return kurtosis_excess() + 3; }
        public abstract double kurtosis_excess();
        public virtual double standard_deviation() { return Math.Sqrt(variance()); }   //either this or variance must be over-ridden to avoid endless loop
        public virtual double hazard(double x)
        {
            double p = cdfc(x);
            double d = pdf(x);
            if (d > p * double.MaxValue) throw new OverflowException();
            if (d == 0) return 0;
            return d / p;
        }
        public virtual double chf(double x) { return -Math.Log(cdfc(x)); }
        protected double inverse_discrete_quantile(double p, double q, double guess, double multiplier, double adder, int max_iter)
        {
            if (p <= pdf(0)) return 0;
            XMath.eps_tolerance tol = new XMath.eps_tolerance(53);
            return do_inverse_discrete_quantile(p, q, guess, multiplier, adder, tol, max_iter);
        }
        double do_inverse_discrete_quantile(double p, double q, double guess, double multiplier, double adder, XMath.eps_tolerance tol, int max_iter)
        {
            distribution_quantile_finder f = new distribution_quantile_finder(this, p, q);
            //
            // Max bounds of the distribution:
            //
            double min_bound, max_bound;
            XMath.pair<double> spt = support();
            min_bound = spt.v1;
            max_bound = spt.v2;

            if (guess > max_bound) guess = max_bound;
            if (guess < min_bound) guess = min_bound;

            double fa = f.next(guess);
            int count = max_iter - 1;
            double fb = fa, a = guess, b = 0.0;

            if (fa == 0) return guess;

            //
            // For small expected results, just use a linear search:
            //
            if (guess < 10)
            {
                b = a;
                while ((a < 10) && (fa * fb >= 0))
                {
                    if (fb <= 0)
                    {
                        a = b;
                        b = a + 1;
                        if (b > max_bound) b = max_bound;
                        fb = f.next(b);
                        --count;
                        if (fb == 0) return b;
                    }
                    else
                    {
                        b = a;
                        a = Math.Max(b - 1, 0.0);
                        if (a < min_bound) a = min_bound;
                        fa = f.next(a);
                        --count;
                        if (fa == 0) return a;
                    }
                }
            }
            //
            // Try and bracket using a couple of additions first, 
            // we're assuming that "guess" is likely to be accurate
            // to the nearest int or so:
            //
            else if (adder != 0)
            {
                //
                // If we're looking for a large result, then bump "adder" up
                // by a bit to increase our chances of bracketing the root:
                //
                //adder = (std::max)(adder, 0.001f * guess);
                if (fa < 0)
                {
                    b = a + adder;
                    if (b > max_bound) b = max_bound;
                }
                else
                {
                    b = Math.Max(a - adder, 0.0);
                    if (b < min_bound) b = min_bound;
                }
                fb = f.next(b);
                --count;
                if (fb == 0) return b;
                if (count > 0 && (fa * fb >= 0))
                {
                    //
                    // We didn't bracket the root, try 
                    // once more:
                    //
                    a = b;
                    fa = fb;
                    if (fa < 0)
                    {
                        b = a + adder;
                        if (b > max_bound) b = max_bound;
                    }
                    else
                    {
                        b = Math.Max(a - adder, 0.0);
                        if (b < min_bound) b = min_bound;
                    }
                    fb = f.next(b);
                    --count;
                }
                if (a > b)
                {
                    XMath.swap(ref a, ref b);
                    XMath.swap(ref fa, ref fb);
                }
            }
            //
            // If the root hasn't been bracketed yet, try again
            // using the multiplier this time:
            //
            if (XMath.sign(fb) == XMath.sign(fa))
            {
                if (fa < 0)
                {
                    //
                    // Zero is to the right of x2, so walk upwards
                    // until we find it:
                    //
                    while (XMath.sign(fb) == XMath.sign(fa))
                    {
                        if (count == 0)
                            throw new Exception(string.Format("Unable to bracket root, last nearest value was {0:G}", b));
                        a = b;
                        fa = fb;
                        b *= multiplier;
                        if (b > max_bound) b = max_bound;
                        fb = f.next(b);
                        --count;
                    }
                }
                else
                {
                    //
                    // Zero is to the left of a, so walk downwards
                    // until we find it:
                    //
                    while (XMath.sign(fb) == XMath.sign(fa))
                    {
                        if (Math.Abs(a) < XMath.min_value)
                        {
                            // Escape route just in case the answer is zero!
                            max_iter -= count;
                            max_iter += 1;
                            return 0;
                        }
                        if (count == 0)
                            throw new Exception(string.Format("Unable to bracket root, last nearest value was {0:G}", a));
                        b = a;
                        fb = fa;
                        a /= multiplier;
                        if (a < min_bound) a = min_bound;
                        fa = f.next(a);
                        --count;
                    }
                }
            }
            max_iter -= count;
            if (fa == 0) return a;
            if (fb == 0) return b;
            //
            // Adjust bounds so that if we're looking for an integer
            // result, then both ends round the same way:
            //
            // adjust_bounds(a, b, tol);    //CGS: Don't think this is needed in our case
            //
            // We don't want zero or denorm lower bounds:
            //
            if (a < XMath.min_value) a = XMath.min_value;
            //
            // Go ahead and find the root:
            //
            XMath.pair<double> r = XMath.toms748_solve(f, a, b, fa, fb, tol, ref count);
            max_iter += count;
            return (r.v1 + r.v2) / 2;
        }
        public virtual double random()
        {
            if (bars_set)
            {
                double translation = (DCL() && DCR()) ? antimode() : 0;
                while (true)
                {
                    double x;
                    int barnum = rand.Next(barcount);
                    bar bar = bars[barnum];
                    if (bar.TailLeft || bar.TailRight)                                                          //All values in the zero block get accepted
                    {
                        x = rand.NextDouble();                                                                  //not an x-value, a cumulative probability instead
                        double plinth = bar.YBottom * bar.Width;                                                //at discontinuities tails are up in the air
                        double zone = bar.Area + bar.TLeft + bar.TRight;                                  
                        x *= (zone + plinth);
                        if((x - plinth) >= bar.TLeft && (x -plinth) <= (bar.TLeft + bar.Area))                  //should only be true for non-discontinuous functions
                            return ((x - bar.TLeft) / bar.Area) * bar.Width + bar.XLeft;                        //Within the rectangle
                        else if (x - plinth < bar.TLeft)
                        {
                            if (translation > 0)
                            {
                                double o = cdf(bar.XLeft);
                                if (x + o > 1) return quantile(x + o - 1);
                                else return quantile(x + o);
                            }
                            return quantile(x);
                        }
                        else
                        {
                            return quantilec(zone + plinth - x);
                        }
                    }
                    else x = bar.XLeft + rand.NextDouble() * bar.Width;
                    if (x >= bar.SafeXLeft && x <= bar.SafeXRight) return translation > 0 && x > 1 ? x - 1 : x;
                    double y = bar.YBottom + rand.NextDouble() * bar.Height;
                    if (translation > 0 && x > 1) x -= 1;
                    if (y <= pdf(x)) return x;
                }
            }
            else
            {
                double p = 0;
                while (p == 0) p = rand.NextDouble();
                double x = p > 0.5 ? quantilec(1 - p) : quantile(p);
                return x;
            }
        }

        bar[] bars;
        protected bool bars_set = false;
        protected bool bars_failed = false;
        class bar
        {
            public double XLeft, XRight, YBottom, YTop, SafeXLeft, SafeXRight;
            public double Width { get { return XRight - XLeft; } }
            public double Height { get { return YTop - YBottom; } }
            public double Area 
            { 
                get 
                {
                    if (double.IsInfinity(Width) && !double.IsInfinity(XLeft) && !double.IsInfinity(XRight))
                        return (-XLeft * Height + XRight * Height);
                    return Height * Width; 
                } 
            }
            public double SafeYTop { get { return SafeXLeft; } set { SafeXLeft = value; } }
            public bool TailLeft, TailRight;
            public double TLeft, TRight;

            #region IBinarySerialize Members

            public void Read(System.IO.BinaryReader r, distribution parent)
            {
                YTop = r.ReadDouble(); ;
                if (parent.LHS())
                {
                    XLeft = r.ReadDouble(); ;
                    SafeXLeft = r.ReadDouble(); ;
                }
                if ((!parent.symmetric() || !parent.LHS()) && parent.RHS())
                {
                    XRight = r.ReadDouble(); ;
                    SafeXRight = r.ReadDouble();
                }
            }

            public void Write(System.IO.BinaryWriter w, distribution parent)
            {
                w.Write(YTop);
                if (parent.LHS())
                {
                    w.Write(XLeft);
                    w.Write(SafeXLeft);
                }
                if ((!parent.symmetric() || !parent.LHS()) && parent.RHS())
                {
                    w.Write(XRight);
                    w.Write(SafeXRight);
                }
            }

            #endregion
        }

        public virtual void setup_bars()
        {
            if(bars_set || bars_failed) return;
            try
            {
                bars = new bar[barcount];
                bars[0] = new bar();
                double A = 1.0 / barcount;                                                              //First guess at A - assumes 100% efficiency so will be too small
                double ymax = max_pdf();
                double ubound = A * 1.5, lbound = A * 0.9;
                double xmode;
                if (unimodal()) xmode = mode();
                else
                {
                    XMath.pair<double> spt = support();
                    if (pdf(spt.v1) > pdf(spt.v2)) xmode = spt.v1;
                    else xmode = spt.v2;
                }

                ziggurat_area_finder zaggy = new ziggurat_area_finder(this, ymax, xmode);
                double fax = zaggy.next(lbound), fbx = zaggy.next(ubound);
                if (Math.Sign(fax) * Math.Sign(fbx) > 0)
                {
                    if (Math.Abs(fax) < 2 * XMath.epsilon || Math.Abs(fbx) < 2 * XMath.epsilon)
                    {
                        bars_set = true;
                        return;
                    }
                    throw new Exception("AAAAAAAAAAAAAAAGH");
                }
                XMath.eps_tolerance tol = new XMath.eps_tolerance(16);
                int max_iter = XMath.max_root_iterations;
                XMath.pair<double> result = XMath.toms748_solve(zaggy, lbound, ubound, fax, fbx, tol, ref max_iter);
                bars_set = true;
              }
            catch (Exception)
            {
                bars_failed = true;
            }
        }

        private void guess_first_bar(double A, bar bar, double ymax, double xmode)
        {
            bar.XLeft = support().v1;
            bar.XRight = support().v2;
            bar.YBottom = 0;
            bar.YTop = ymax / barcount;                                                             //First guess at height of base bar
            if (!LHS())                                                                             //LHS is vertical 
            {
                bar.SafeXLeft = bar.XLeft;
                bar.TailLeft = false;
                bar.TLeft = 0;
            }
            else if (!tail_left())                                                                  //LHS starts at a finite value - no tail to worry about
            {
                bar.TailLeft = false;
                bar.TLeft = 0;
                if (bar.YTop <= pdf(bar.XLeft)) bar.SafeXLeft = bar.XLeft;                          //Still vertical at this point - curve starts higher up y axis
                else bar.SafeXLeft = pdf_inv(bar.YTop, false);
            }
            else                                                                                    //Tail should be put outside the rectangle
            {
                bar.TailLeft = true;
            }
            if (!RHS())                                                                             //RHS is vertical 
            {
                bar.SafeXRight = bar.XRight;
                bar.TailRight = false;
                bar.TRight = 0;
            }
            else if (!tail_right())                                                                   //RHS ends at a finite value - no tail to worry about
            {
                bar.TailRight = false;
                bar.TRight = 0;
                if (bar.YTop <= pdf(bar.XRight)) bar.SafeXRight = bar.XRight;
                else bar.SafeXRight = pdf_inv(bar.YTop, true);
            }
            else                                                                                    //Tail should be put outside the rectangle
            {
                bar.TailRight = true;
            }

            ziggurat_bar_finder ziggy = new ziggurat_bar_finder(this, bar, A, xmode);
            double lbound = double.Epsilon, ubound = bar.YTop;
            double fax = ziggy.next(lbound), fbx = ziggy.next(ubound);
            if (Math.Sign(fax) * Math.Sign(fbx) > 0)
            {
                if (Math.Abs(fax) < 2 * XMath.epsilon || Math.Abs(fbx) < 2 * XMath.epsilon)
                {
                    bars_set = true;
                    return;
                }
                throw new Exception("AAAAAAAAAAAAAAAGH");
            }
            XMath.eps_tolerance tol = new XMath.eps_tolerance(16);
            int max_iter = XMath.max_root_iterations;
            XMath.pair<double> result = XMath.toms748_solve(ziggy, lbound, ubound, fax, fbx, tol, ref max_iter);
        }

        class ziggurat_bar_finder : XMath.iterand<double, double>
        {
            distribution dist;
            bar bar;
            double target_area, xmode;

            public ziggurat_bar_finder(distribution d, bar b, double a, double m)
            {
                dist = d;
                bar = b;
                target_area = a;
                xmode = m;
            }

            public override double next(double ytop)
            {
                bar.YTop = ytop;
                if (bar.TailLeft)
                {
                    bar.XLeft = bar.SafeXLeft = dist.pdf_inv(bar.YTop, false);
                    bar.TLeft = dist.cdf(bar.SafeXLeft);
                    //bar.YTop = dist.pdf(bar.SafeXLeft);                                                  //Apply a correction in case pdf_inv is not accurate
                }
                if (dist.symmetric())
                {
                    bar.XRight = xmode + (xmode - bar.XLeft);
                    if (double.IsInfinity(bar.XRight)) bar.XRight = double.MaxValue;
                    bar.SafeXRight = xmode + (xmode - bar.SafeXLeft);
                    if (double.IsInfinity(bar.SafeXRight)) bar.SafeXRight = double.MaxValue;
                    bar.TRight = bar.TLeft;
                }
                else if (bar.TailRight)
                {
                    bar.XRight = bar.SafeXRight = dist.pdf_inv(bar.YTop, true);
                    bar.TRight = dist.cdfc(bar.SafeXRight);
                }
                double A0 = bar.TLeft + bar.TRight + bar.Area;
                return target_area - A0;
            }
        }

        class ziggurat_area_finder : XMath.iterand<double, double>
        {
            distribution dist;
            double target;
            double xmode;
            double translation;

            public ziggurat_area_finder(distribution d, double ymax, double xm)
            {
                dist = d;
                target = ymax;
                xmode = xm;
                translation = 0;                                                        //Amount to move x through to turn U=shaped beta distribution into volcano-shaped one.
            }

            public override double next(double A)
            {
                int bartop = dist.barcount, barbot = 1;
                if (dist.DCL() || dist.DCR())
                {
                    if(dist.bars[--bartop] == null) dist.bars[bartop] = new bar();
                    dist.find_ymax(A, out target, out translation, dist.bars[bartop]);
                }
                int i;
                if (dist.DCL() && dist.DCR())
                {
                    int rtop = (int)(dist.pdf(dist.antimode()) / A);
                    for (i = 0; i < rtop; i++)
                    {
                        if (dist.bars[i] == null) dist.bars[i] = new bar();
                        dist.bars[i].SafeXLeft = dist.bars[i].XLeft = translation;
                        dist.bars[i].SafeXRight = dist.bars[i].XRight = translation + 1;
                        dist.bars[i].YBottom = i == 0 ? 0 : dist.bars[i - 1].YTop;
                        dist.bars[i].YTop = A + dist.bars[i].YBottom;
                    }
                    barbot = i;
                }
                else dist.guess_first_bar(A, dist.bars[0], target, xmode);
                for (i = barbot; i < bartop; i++)
                {
                    if (dist.bars[i] == null) dist.bars[i] = new bar();
                    dist.bars[i].XLeft = dist.bars[i - 1].SafeXLeft;
                    dist.bars[i].XRight = dist.bars[i - 1].SafeXRight;
                    dist.bars[i].YBottom = dist.bars[i - 1].YTop;
                    dist.bars[i].YTop = A / dist.bars[i].Width + dist.bars[i].YBottom;
                    if (dist.bars[i].YTop > target) break;

                    if (dist.DCL() && dist.DCR())
                    {
                        dist.bars[i].SafeXLeft = dist.pdf_inv(dist.bars[i].YTop, true);
                        dist.bars[i].SafeXRight = dist.pdf_inv(dist.bars[i].YTop, false) + 1;
                    }
                    else
                    {
                        if (dist.LHS() && dist.bars[i].YTop > dist.pdf(dist.bars[0].XLeft)) dist.bars[i].SafeXLeft = dist.pdf_inv(dist.bars[i].YTop, false);
                        else dist.bars[i].SafeXLeft = dist.bars[i].XLeft;
                        if (dist.symmetric()) dist.bars[i].SafeXRight = xmode + (xmode - dist.bars[i].SafeXLeft);
                        else if (dist.RHS() && dist.bars[i].YTop > dist.pdf(dist.bars[0].XRight)) dist.bars[i].SafeXRight = dist.pdf_inv(dist.bars[i].YTop, true);
                        else dist.bars[i].SafeXRight = dist.bars[i].XRight;
                    }
                }
                if (i < bartop) return ((double)i / (double)bartop - 1.0);
                return target - dist.bars[i - 1].YTop;
            }
        }

        class distribution_invpdf_finder : XMath.iterand<double, double>
        {
            distribution dist;
            double target;

            public distribution_invpdf_finder(distribution d, double p)
            {
                dist = d;
                target = p;
            }

            public override double next(double x)
            {
                return target - dist.pdf(x);
            }
        }

        class distribution_quantile_finder : XMath.iterand<double, double>
        {
            distribution dist;
            double target;
            bool comp;

            public distribution_quantile_finder(distribution d, double p, double q)
            {
                dist = d;
                target = p < q ? p : q;
                comp = p < q ? false : true;
            }

            public override double next(double x)
            {
                return comp ? target - dist.cdfc(x) : dist.cdf(x) - target;
            }
        }

        private double leap_left(double x, double limit)
        {
            double result;
            if (x > Math.E) result = Math.Exp(Math.Log(x) * 0.5);
            else if (x > 1)
            {
                if (limit >= 1) result = Math.Exp(Math.Log(x) * 0.5); //Keep shrinking - close to other bound
                else result = 1;
            }
            else if (x > 0)
            {
                if (limit >= 0 && x > double.Epsilon)
                {
                    if (x == 1) result = 1 - XMath.epsilon;
                    else result = Math.Exp(Math.Log(x) * 2);             //Shrink some more - exponent is now negative so double it to shrink further
                }
                else result = 0;
            }
            else if (x > -1)
            {
                if (limit >= -1)
                {
                    if (x == 0) result = -double.Epsilon;
                    else result = -Math.Exp(Math.Log(-x) * 0.5);             //Grow now in negative direction
                }
                else result = -1;
            }
            else if (x > -Math.E)
            {
                if (limit >= -Math.E)
                {
                    if (x == -1) result = -(1 + XMath.epsilon);
                    else result = -Math.Exp(Math.Log(-x) * 2);             //Grow now in negative direction
                }
                else result = -Math.E;
            }
            else result = -Math.Exp(Math.Log(-x) * 2);
            return result;
        }

        private double leap_right(double x, double limit)
        {
            double result;
            if (x < -Math.E) result = -Math.Exp(Math.Log(-x) * 0.5);
            else if (x < -1)
            {
                if (limit <= -1) result = -Math.Exp(Math.Log(-x) * 0.5);
                else  result = -1;
            }
            else if (x < 0)
            {
                if (limit <= 0)
                {
                    if (x == -1) result = -1 + XMath.epsilon;
                    else result = -Math.Exp(Math.Log(-x) * 2);
                }
                else result = 0;
            }
            else if (x < 1)
            {
                if (limit <= 1)
                {
                    if (x == 0) result = XMath.epsilon;
                    else result = Math.Exp(Math.Log(x) * 0.5);
                }
                else result = 1;
            }
            else if (x < Math.E)
            {
                if (limit <= Math.E)
                {
                    if (x == 1) result = 1 + XMath.epsilon;
                    else result = Math.Exp(Math.Log(x) * 2);
                }
                else result = Math.E;
            }
            else result = Math.Exp(Math.Log(x) * 2);
            return result;
        }

        protected double find_pdf_inv(double p, double lbound, double ubound, bool increasing)
        {

            ///Try to narrow down the bounds so the iteration has some chance
            ///Principle is to move them exponentially towards each other until they cross the root
            ///When this happens, the 'overstep' may be a suitable value for the opposite bound.
            ///Assumption at this point is that there are no maxima or minima between the two bounds
            ///But we need to be careful when 0, 1 or e lie in the interesting region

            if (ubound == lbound) return ubound;
            if (ubound < lbound) throw new Exception("Inverse PDF internal exception: upper bound and lower bound are the wrong way round");
            if (increasing && pdf(lbound) > pdf(ubound)) throw new Exception("Inverse PDF internal exception: sequence specified as increasing but appears to be decreasing");
            else if (!increasing && pdf(lbound) < pdf(ubound)) throw new Exception("Inverse PDF internal exception: sequence specified as decreasing but appears to be increasing");

            double min = lbound, max = ubound;  //Can't go outside the original boundaries
            double l, u;

            do
            {
                l = leap_right(lbound, ubound);
                if (l > max) break;
                if ((increasing && pdf(l) < p) || (!increasing && pdf(l) > p)) lbound = l;
                else if (ubound > l) ubound = l;
            } while (lbound == l);

            do
            {
                u = leap_left(ubound, lbound);
                if (u < min) break;
                if ((increasing && pdf(u) > p) || (!increasing && pdf(u) < p)) ubound = u;
                else if (lbound < u) lbound = u;
            } while (ubound == u);

            //Might already have done enough
            if (ubound - lbound < 2 * XMath.epsilon) return lbound + 0.5 * (ubound - lbound);

            distribution_invpdf_finder finder = new distribution_invpdf_finder(this, p);
            double fax = finder.next(lbound), fbx = finder.next(ubound);
            if (Math.Sign(fax) * Math.Sign(fbx) > 0)
            {
                if (Math.Abs(fax - fbx) < 2 * XMath.epsilon) return lbound + 0.5 * (ubound - lbound);
                if (Math.Abs(fax) < 2 * XMath.epsilon) return lbound;
                if (Math.Abs(fbx) < 2 * XMath.epsilon) return ubound;
                //return 0;
                throw new Exception("AAAAAAAAAAAAAAAGH");
            }
            XMath.eps_tolerance tol = new XMath.eps_tolerance(53);
            int max_iter = XMath.max_root_iterations;
            XMath.pair<double> result = XMath.toms748_solve(finder, lbound, ubound, fax, fbx, tol, ref max_iter);
            if (max_iter < XMath.max_root_iterations || (Math.Abs(result.v1 - result.v2) < 2 * double.Epsilon)) 
                return result.v1 + 0.5 * (result.v2 - result.v1);
            else
                result = XMath.bisect(finder, lbound, ubound, tol, ref max_iter);
            if (max_iter < XMath.max_root_iterations || (Math.Abs(result.v1 - result.v2) < 2 * double.Epsilon))
                return result.v1 + 0.5 * (result.v2 - result.v1);
            else
                throw new Exception("Inverse PDF: could not get series to converge");
        }

        internal double vertical_tail_area(double y)
        {
            double x, total, rect, result = 0;

            if (y == 0) return 1;
            if (DCL() && DCR())
            {
                x = pdf_inv(y, false);                          //U-shaped beta - RHS and LHS get swapped around
                total = cdf(x);
                rect = (x - support().v1) * pdf(x);
                result += total - rect;
                x = pdf_inv(y, true);
                total = cdfc(x);
                rect = (support().v2 - x) * pdf(x);
                result += total - rect;
            }
            else if (DCL())
            {
                x = pdf_inv(y, true);
                total = cdf(x);
                rect = (x - support().v1) * y;
                result += total - rect;
            }
            else if (DCR())
            {
                x = pdf_inv(y, false);
                total = cdfc(x);
                rect = (support().v2 - x) * pdf(x);
                result += total - rect;
            }
            return result;
        }

        class ymax_finder : XMath.iterand<double, double>
        {
            double area;
            distribution dist;

            public ymax_finder(distribution d, double A)
            {
                area = A;
                dist = d;
            }

            public override double next(double y)
            {
                return(area - dist.vertical_tail_area(y));
            }
        }
        
        private void find_ymax(double area, out double target, out double translation, bar bar)
        {
            double lbound, ubound, fax, fbx, ymax = 0;

            translation = 0;
            if (DCL() && DCR())
            {
                translation = antimode();
                ubound = Math.Max(pdf(quantile(area / 2)), pdf(quantilec(area / 2)));
            }
            else if (DCL()) ubound = pdf(quantile(area));
            else ubound = pdf(quantilec(area));
            if (DCL() && DCR()) lbound = pdf(antimode());                                       //U -shaped beta
            else lbound = double.Epsilon;
            ymax_finder ymf = new ymax_finder(this, area);
            fax = ymf.next(lbound);
            fbx = ymf.next(ubound);
            if (Math.Sign(fax) * Math.Sign(fbx) > 0)
            {
                if (Math.Abs(fax) <= 2 * XMath.epsilon) ymax = lbound;
                else if (Math.Abs(fbx) <= 2 * XMath.epsilon) ymax = ubound;
                else throw new Exception("AAAAAAAAAAAAAAGH");
            }
            else
            {
                XMath.eps_tolerance tol = new XMath.eps_tolerance(53);
                int max_iter = XMath.max_root_iterations;
                XMath.pair<double> bounds = XMath.toms748_solve(ymf, lbound, ubound, fax, fbx, tol, ref max_iter);
                ymax = (bounds.v1 + bounds.v2) / 2;
            }
            target = ymax;
            bar.YTop = bar.YBottom = ymax;
            if (DCL() && DCR())
            {
                bar.XLeft = pdf_inv(ymax, true);          //implement wrap-around
                bar.XRight = pdf_inv(ymax, false) + 1;        //implement wrap-around
                bar.TailLeft = true;
                bar.TLeft = vertical_tail_area(ymax);
            }
            else if (DCL())
            {
                bar.XLeft = support().v1;
                bar.XRight = pdf_inv(ymax, true);
                bar.TailLeft = true;
                bar.TLeft = vertical_tail_area(ymax);
            }
            else
            {
                bar.XLeft = pdf_inv(ymax, false);
                bar.XRight = support().v2;
                bar.TailRight = true;
                bar.TRight = vertical_tail_area(ymax);
            }


        }

        class generic_quantile_finder : XMath.iterand<double, double>
        {
            distribution dist;
            double target;
            bool comp;

            public generic_quantile_finder(distribution d, double t, bool c)
            {
                dist = d;
                target = t;
                comp = c;
            }

            public override double next(double x)
            {
                return comp ? target - dist.cdfc(x) : dist.cdf(x) - target;
            }

        }

        protected double generic_quantile(double p, double guess, bool comp)
        {
            if (p == 0) return comp ? range().v2 : range().v1;
            if (p == 1) return comp ? range().v1 : range().v2;

            generic_quantile_finder f = new generic_quantile_finder(this, p, comp);
            XMath.eps_tolerance tol = new XMath.eps_tolerance(50);
            int max_iter = XMath.max_root_iterations;
            XMath.pair<double> ir = XMath.bracket_and_solve_root(f, guess, 2.0, true, tol, ref max_iter);
            double result = ir.v1 + (ir.v2 - ir.v1) / 2;
            if (max_iter >= XMath.max_root_iterations) throw new Exception(this.GetType().Name + ": quantile: unable to locate solution in a reasonable time");
            return result;
        }

        class pdf_maximizer : XMath.iterand<double, double>
        {
            distribution dist;

            public pdf_maximizer(distribution d)
            {
                dist = d;
            }

            public override double next(double x)
            {
                return -dist.pdf(x);
            }
        }

        class pdf_minimizer : XMath.iterand<double, double>
        {
            distribution dist;

            public pdf_minimizer(distribution d)
            {
                dist = d;
            }

            public override double next(double x)
            {
                return dist.pdf(x);
            }
        }

        protected double generic_find_mode(double guess, double step)
        {
            double maxval;
            double upper_bound = guess;
            double lower_bound;
            double v = pdf(guess);
            if (v == 0) throw new Exception(string.Format("{0:G}: Mode: Could not locate a starting location for the search for the mode, original guess was {1:G}", this.GetType().Name, guess));
            do
            {
                maxval = v;
                if (step != 0) upper_bound += step;
                else upper_bound *= 2;
                v = pdf(upper_bound);
            }
            while (maxval < v);

            lower_bound = upper_bound;
            do
            {
                maxval = v;
                if (step != 0)
                    lower_bound -= step;
                else
                    lower_bound /= 2;
                v = pdf(lower_bound);
            }
            while (maxval < v);

            int max_iter = XMath.max_root_iterations;
            double result = brent_find_minima(new pdf_maximizer(this), lower_bound, upper_bound, ref max_iter).v1;
            if (max_iter >= XMath.max_root_iterations)
                throw new Exception(string.Format("{0:G}: Mode: Unable to locate solution in a reasonable time, current best guess is {1:G}", this.GetType().Name, guess));
            return result;
        }

        XMath.pair<double> brent_find_minima(XMath.iterand<double, double> f, double min, double max, ref int max_iter)
        {
            int bits = 26;
            double tolerance = XMath.ldexp(1.0, 1 - bits);
            double x;  // minima so far
            double w;  // second best point
            double v;  // previous value of w
            double u;  // most recent evaluation point
            double delta;  // The distance moved in the last step
            double delta2; // The distance moved in the step before last
            double fu, fv, fw, fx;  // function evaluations at u, v, w, x
            double mid; // midpoint of min and max
            double fract1, fract2;  // minimal relative movement in x

            double golden = 0.3819660;  // golden ratio, don't need too much precision here!

            x = w = v = max;
            fw = fv = fx = f.next(x);
            delta2 = delta = 0;

            int count = max_iter;

            do
            {
                // get midpoint
                mid = (min + max) / 2;
                // work out if we're done already:
                fract1 = tolerance * Math.Abs(x) + tolerance / 4;
                fract2 = 2 * fract1;
                if (Math.Abs(x - mid) <= (fract2 - (max - min) / 2))
                    break;

                if (Math.Abs(delta2) > fract1)
                {
                    // try and construct a parabolic fit:
                    double r = (x - w) * (fx - fv);
                    double q = (x - v) * (fx - fw);
                    double p = (x - v) * q - (x - w) * r;
                    q = 2 * (q - r);
                    if (q > 0)
                        p = -p;
                    q = Math.Abs(q);
                    double td = delta2;
                    delta2 = delta;
                    // determine whether a parabolic step is acceptible or not:
                    if ((Math.Abs(p) >= Math.Abs(q * td / 2)) || (p <= q * (min - x)) || (p >= q * (max - x)))
                    {
                        // nope, try golden section instead
                        delta2 = (x >= mid) ? min - x : max - x;
                        delta = golden * delta2;
                    }
                    else
                    {
                        // whew, parabolic fit:
                        delta = p / q;
                        u = x + delta;
                        if (((u - min) < fract2) || ((max - u) < fract2))
                            delta = (mid - x) < 0 ? (double)-Math.Abs(fract1) : (double)Math.Abs(fract1);
                    }
                }
                else
                {
                    // golden section:
                    delta2 = (x >= mid) ? min - x : max - x;
                    delta = golden * delta2;
                }
                // update current position:
                u = (Math.Abs(delta) >= fract1) ? x + delta : (delta > 0 ? x + Math.Abs(fract1) : x - Math.Abs(fract1));
                fu = f.next(u);
                if (fu <= fx)
                {
                    // good new point is an improvement!
                    // update brackets:
                    if (u >= x)
                        min = x;
                    else
                        max = x;
                    // update control points:
                    v = w;
                    w = x;
                    x = u;
                    fv = fw;
                    fw = fx;
                    fx = fu;
                }
                else
                {
                    // Oh dear, point u is worse than what we have already,
                    // even so it *must* be better than one of our endpoints:
                    if (u < x)
                        min = u;
                    else
                        max = u;
                    if ((fu <= fw) || (w == x))
                    {
                        // however it is at least second best:
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    }
                    else if ((fu <= fv) || (v == x) || (v == w))
                    {
                        // third best:
                        v = u;
                        fv = fu;
                    }
                }

            } while (--count > 0);

            max_iter -= count;

            return new XMath.pair<double>(x, fx);
        }

        protected double generic_find_mode_01(double guess)
        {
            // Need to begin by bracketing the maxima of the PDF:
            //
            double maxval;
            double upper_bound = guess;
            double lower_bound;
            double v = pdf(guess);
            do
            {
                maxval = v;
                upper_bound = 1 - (1 - upper_bound) / 2;
                if (upper_bound == 1) return 1;
                v = pdf(upper_bound);
            } while (maxval < v);

            lower_bound = upper_bound;
            do
            {
                maxval = v;
                lower_bound /= 2;
                if (lower_bound < double.Epsilon) return 0;
                v = pdf(lower_bound);
            }
            while (maxval < v);

            int max_iter = XMath.max_root_iterations;

            double result = brent_find_minima(new pdf_maximizer(this), lower_bound, upper_bound, ref max_iter).v1;
            if (max_iter >= XMath.max_root_iterations)
                throw new Exception(string.Format("{0:G}: Mode: Unable to locate solution in a reasonable time, current best guess is {1:G}", this.GetType().Name, guess));
            return result;
        }
 
        protected double generic_find_antimode_01(double guess)
        {
            // Need to begin by bracketing the minima of the PDF:
            //
            double minval;
            double upper_bound = guess;
            double lower_bound;
            double v = pdf(guess);
            do
            {
                minval = v;
                upper_bound = 1 - (1 - upper_bound) / 2;
                if (upper_bound == 1) return 1;
                v = pdf(upper_bound);
            } while (minval > v);

            lower_bound = upper_bound;
            do
            {
                minval = v;
                lower_bound /= 2;
                if (lower_bound < double.Epsilon) return 0;
                v = pdf(lower_bound);
            }
            while (minval > v);

            int max_iter = XMath.max_root_iterations;

            double result = brent_find_minima(new pdf_minimizer(this), lower_bound, upper_bound, ref max_iter).v1;
            if (max_iter >= XMath.max_root_iterations)
                throw new Exception(string.Format("{0:G}: Mode: Unable to locate solution in a reasonable time, current best guess is {1:G}", this.GetType().Name, guess));
            return result;
        }

        #region IBinarySerialize Members

        public void Read(System.IO.BinaryReader r)
        {
            bars_set = r.ReadBoolean();
            if (bars_set)
            {
                bars[0].TailLeft = r.ReadBoolean();
                if (bars[0].TailLeft) bars[0].TLeft = r.ReadDouble();
                bars[0].TailRight = r.ReadBoolean();
                if (bars[0].TailRight) bars[0].TRight = r.ReadDouble();
                if (!LHS()) bars[0].XLeft = r.ReadDouble();
                if (!RHS()) bars[0].XRight = r.ReadDouble();
                bars[0].YBottom = 0;
                for (int i = 0; i < barcount; i++)
                {
                    bars[i].Read(r, this);
                    if (i > 0)
                    {
                        if (!LHS()) bars[i].XLeft = bars[i].SafeXLeft = bars[i - 1].XLeft;
                        if (!RHS()) bars[i].XRight = bars[i].SafeXRight = bars[i - 1].XRight;
                        bars[i].YBottom = bars[i - 1].YTop;
                    }
                }
            }
        }

        public void Write(System.IO.BinaryWriter w)
        {
            w.Write(bars_set);
            if (bars_set)
            {
                w.Write(bars[0].TailLeft);
                if (bars[0].TailLeft) w.Write(bars[0].TLeft);
                w.Write(bars[0].TailRight);
                if (bars[0].TailRight) w.Write(bars[0].TRight);
                if (!LHS()) w.Write(bars[0].XLeft);
                if (!RHS()) w.Write(bars[0].XRight);
                for (int i = 0; i < barcount; i++) bars[i].Write(w, this);
            }
        }

        #endregion

    }



}
