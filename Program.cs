
using System;
using System.Collections.Generic;
using System.Linq;
namespace Optimization
{
    class CIndividual
    {
        public double Fitness = 0;
        public double[] Points ;
    }
    class CHE
    {

        double NOFE;
        private int var;
        private bool SUCCES = false;
        private readonly double GlobOptimization;
        private readonly double stop_Condition;
        private int iteration;
        private double upperLimit;
        private double lowerLimit;
        private Random generater ;
        private CIndividual[] mPopulation;
        private int n = 0;
        private CIndividual Best_Sol;
        double Radius = 0;
        List<double> LastGeneration;
        public CHE(int var,
                            int n,
                            int iteration,
                            double upperLimit,
                            double lowerLimit,
                            double globalOptimization,
                            double stopCondition
                            )
        {

            this.n = n;
            this.var = var;
            this.iteration = iteration;
            this.upperLimit = upperLimit;
            this.lowerLimit = lowerLimit;
            this.GlobOptimization = globalOptimization;
            this.stop_Condition = stopCondition;
            this.SUCCES = false;
            this.mPopulation = new CIndividual[this.n];
            this.generater = new Random();
            this.Radius = (this.upperLimit - this.lowerLimit) / 1800; // 1800 soucified by user, could be any number.
            NOFE = 0;
            LastGeneration = new List<double>();

        }

        public void run()
        {
            Best_Sol = new CIndividual();
            Best_Sol.Points = new double[var];
            Best_Sol.Fitness = double.MaxValue;

            for (int i = 0; i < this.n; i++)
            {
                mPopulation[i] = new CIndividual();
                mPopulation[i].Points = new double[this.var];

                for (int j = 0; j < this.var; j++)
                {
                    mPopulation[i].Points[j] = generater.NextDouble() * (this.upperLimit - this.lowerLimit) + this.lowerLimit;
                }

                mPopulation[i].Fitness = calculateFitness(mPopulation[i].Points);

                if (mPopulation[i].Fitness < Best_Sol.Fitness)
                {
                    Best_Sol.Fitness = mPopulation[i].Fitness;
                    mPopulation[i].Points.CopyTo(Best_Sol.Points, 0);
                }
            }

            this.Succes = false;
            int itr = 1;
            double[] x_elite = new double[this.var];
            double bestAv = Best_Sol.Points.Average();
            for (itr = 0; itr < iteration; itr++)//main loop
            {
                double result = 0;
                if (Best_Sol.Fitness >= this.GlobOptimization)
                    result = Best_Sol.Fitness - this.GlobOptimization;
                else
                    result = this.GlobOptimization - Best_Sol.Fitness;

                if (result <= this.stop_Condition)
                {
                    this.Succes = true;
                    break;
                }

                for (int i = 0; i < this.n / 2; i++)
                {

                    for (int j = 0; j < this.var; j++)
                    {
                        double x = (bestAv + mPopulation[i].Points[j]) / 2;
                        x_elite[j] = (generater.NextDouble() * 2 - 1) * (x + Best_Sol.Points[j]);

                        if (x_elite[j] > this.upperLimit) x_elite[j] = this.upperLimit;
                        if (x_elite[j] < this.lowerLimit) x_elite[j] = this.lowerLimit;
                    }
                    double f = calculateFitness(x_elite);

                    if (f < Best_Sol.Fitness)
                    {
                        Best_Sol.Fitness = f;
                        x_elite.CopyTo(Best_Sol.Points, 0);
                        bestAv = Best_Sol.Points.Average();
                    }
                    if (f < mPopulation[i].Fitness)
                    {
                        mPopulation[i].Fitness = f;
                        x_elite.CopyTo(mPopulation[i].Points, 0);
                    }

                }

                for (int i = this.n / 2; i < this.n; i++)
                {
                    for (int j = 0; j < this.var; j++)
                    {
                        if (generater.NextDouble() >= 0.5)
                            x_elite[j] = (generater.NextDouble()) * (mPopulation[i].Points[j] - Best_Sol.Points.Min()) + Best_Sol.Points.Min();
                        else
                            x_elite[j] = (generater.NextDouble()) * (mPopulation[i].Points[j] - Best_Sol.Points[j]) + Best_Sol.Points[j];

                        if (x_elite[j] > this.upperLimit) x_elite[j] = this.upperLimit;
                        if (x_elite[j] < this.lowerLimit) x_elite[j] = this.lowerLimit;
                    }

                    double f = calculateFitness(x_elite);

                    if (f < Best_Sol.Fitness)
                    {
                        Best_Sol.Fitness = f;
                        x_elite.CopyTo(Best_Sol.Points, 0);
                        bestAv = Best_Sol.Points.Average();
                    }
                    if (f < mPopulation[i].Fitness)
                    {
                        mPopulation[i].Fitness = f;
                        x_elite.CopyTo(mPopulation[i].Points, 0);
                    }

                    Best_Sol.Points.CopyTo(x_elite, 0);
                    int loc = generater.Next(0, this.var);

                    if (generater.NextDouble() <= 0.5)
                        x_elite[loc] = generater.NextDouble() * ((Best_Sol.Points[loc] + Radius) - (Best_Sol.Points[loc] - Radius)) + (Best_Sol.Points[loc] - Radius);
                    else
                        x_elite[loc] = generater.NextDouble() * (this.upperLimit - this.lowerLimit) + this.lowerLimit;

                    if (x_elite[loc] > this.upperLimit) x_elite[loc] = this.upperLimit;
                    if (x_elite[loc] < this.lowerLimit) x_elite[loc] = this.lowerLimit;

                    f = calculateFitness(x_elite);

                    if (f < Best_Sol.Fitness)
                    {
                        Best_Sol.Fitness = f;
                        x_elite.CopyTo(Best_Sol.Points, 0);
                        bestAv = Best_Sol.Points.Average();
                    }
                }

            }//end main iteration

            for (int i = 0; i < this.n; i++)
                LastGeneration.Add(mPopulation[i].Fitness);

        }
        public double GetNOFE()
        {
            return NOFE;
        }

        public double BestFintness()
        {
            return Best_Sol.Fitness;
        }

        public double[] GeBestPoints()
        {

            return Best_Sol.Points;

        }
        public bool Succes
        {
            get { return SUCCES; }
            set { SUCCES = value; }
        }
        public double[] GetLastGeneration()
        {
            double[] temp = new double[LastGeneration.Count];
            temp = LastGeneration.ToArray();
            return temp;
        }

        private double calculateFitness(double[] XD)
        {
            NOFE++;
            double fit = double.MaxValue;
            fit = Function1(XD);
            //switch (FunctionIndex)
            //{
            //    case 0: fit = Function1(XD); break;
            //    case 1: fit = Function2(XD); break;
            //    case 2: fit = Function3(XD); break;
            //    case 3: fit = Function4(XD); break;
            //    case 4: fit = Function5(XD); break;
            //    case 5: fit = Function6(XD); break;
            //    case 6: fit = Function7(XD); break;
            //    case 7: fit = Function8(XD); break;
            //    case 8: fit = Function9(XD); break;
            //    case 9: fit = Function10(XD); break;
            //    case 10: fit = Function11(XD); break;
            //    case 11: fit = Function12(XD); break;
            //    case 12: fit = Function13(XD); break;
            //    case 13: fit = Function14(XD); break;
            //    case 14: fit = Function15(XD); break;
            //    case 15: fit = Function16(XD); break;
            //    case 16: fit = Function17(XD); break;
            //    case 17: fit = Function18(XD); break;
            //    case 18: fit = Function19(XD); break;
            //    case 19: fit = Function20(XD); break;
            //    case 20: fit = Function21(XD); break;
            //    case 21: fit = Function22(XD); break;
            //    case 22: fit = Function23(XD); break;
            //    case 23: fit = Function24(XD); break;
            //    case 24: fit = Function25(XD); break;
            //    case 25: fit = Function26(XD); break;
            //    case 26: fit = Function27(XD); break;
            //    case 27: fit = Function28(XD); break;
            //    case 28: fit = Function29(XD); break;
            //    case 29: fit = Function30(XD); break;
            //    case 30: fit = Function31(XD); break;
            //    case 31: fit = Function32(XD); break;
            //    case 32: fit = Function33(XD); break;
            //    case 33: fit = Function34(XD); break;
            //    case 34: fit = Function35(XD); break;
            //    case 35: fit = Function36(XD); break;
            //    case 36: fit = Function37(XD); break;
            //    case 37: fit = Function38(XD); break;
            //    case 38: fit = Function39(XD); break;
            //    case 39: fit = Function40(XD); break;
            //    case 40: fit = Function41(XD); break;
            //    case 41: fit = Function42(XD); break;
            //    case 42: fit = Function43(XD); break;
            //    case 43: fit = Function44(XD); break;
            //    case 44: fit = Function45(XD); break;
            //    case 45: fit = Function46(XD); break;
            //    case 46: fit = Function47(XD); break;
            //    case 47: fit = Function48(XD); break;
            //    case 48: fit = Function49(XD); break;
            //    case 49: fit = Function50(XD); break;

            //}
            return fit;
        }// end funct


        /****************************************************
          * Sphere 50 D,  range = [-100, 100] 
          *   Min = 0,     x = 0,0....0
          *    
          * ***************************************************/
        private double Function1(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
                sum += Math.Pow(XD[i], 2);
            return sum;
        }

        // Quartic Noise 50 D, range -1.28, 1.28 min = 0
        private double Function2(double[] XD)
        {

            //Ok
            double sum = 0.0;
            for (int i = 0; i < XD.Length; i++)
            {

                sum += (((i + 1) * Math.Pow(XD[i], 4)));
            }

            return sum + generater.NextDouble();
        }
        //Powell Sum 50 D, range [-1, 1] min=0
        private double Function3(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            int n = 2;
            for (int i = 0; i < XD.Length; i++)
            {
                double abs = Math.Abs(XD[i]);
                sum += (Math.Pow(abs, n));
                n++;
            }
            return sum;
        }

        //Schwefel 2.20 D50  range -100,100 min=0
        private double Function4(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                double abs = Math.Abs(XD[i]);
                sum += abs;
            }
            return sum;
        }

        //Schwefel 2.21 D50  range -100,100 min=0
        private double Function5(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double[] temp = new double[XD.Length];

            for (int i = 0; i < XD.Length; i++)
                temp[i] = Math.Abs(XD[i]);
            sum = temp.Max();
            return sum;
        }
        //Step 50 D , range [-100, 100] min = 0
        private double Function6(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                sum += Math.Pow(Math.Floor(XD[i] + 0.5), 2);
            }
            return sum;
        }



        /****************************************************
   * stepint 50 D
   * range = [-5.12,  5.12]   Min = 25-6n,     
   * 
   * 
   * ***************************************************/
        private double Function7(double[] XD)
        {
            //Ok

            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
                sum += Math.Floor(XD[i]);
            sum += 25;
            return sum;
        }

        //Schwefel 1.2 D 50  -100, 100  min = 0
        private double Function8(double[] XD)
        {
            double sum = 0;
            double p1;

            for (int i = 0; i < XD.Length; i++)
            {
                p1 = 0;
                for (int j = 0; j <= i; j++)
                    p1 += XD[j];
                sum += Math.Pow(p1, 2);

            }
            return sum;
        }

        // Schwefel 2.22 D 50 -100, 100
        private double Function9(double[] XD)
        {
            double sum = 0;
            double p1 = 0;
            double p2 = 1;

            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Abs(XD[i]);
                p2 *= Math.Abs(XD[i]);

            }
            sum = p1 + p2;
            return sum;
        }

        // Schwefel 2.23 D 50 -10, 10
        private double Function10(double[] XD)
        {
            double sum = 0;

            for (int i = 0; i < XD.Length; i++)
            {
                sum += Math.Pow(XD[i], 10);
            }
            return sum;
        }
        //Rosenbrock 50 D,  -30,30  fit=0,  X=1,1,1...1
        private double Function11(double[] XD)
        {
            double sum = 0;
            for (int i = 0; i < XD.Length - 1; i++)
            {
                double p1 = XD[i + 1] - Math.Pow(XD[i], 2);
                double p2 = Math.Pow(XD[i] - 1, 2);
                sum += 100 * Math.Pow(p1, 2) + p2;

            }
            return sum;
        }

        //Brown 50 D,  -1,4  fit=0,  
        private double Function12(double[] XD)
        {
            double sum = 0;
            for (int i = 0; i < XD.Length - 1; i++)
            {
                double p1 = Math.Pow(XD[i], 2);
                double p2 = Math.Pow(XD[i + 1], 2);
                sum += Math.Pow(p1, p2 + 1) + Math.Pow(p2, p1 + 1);

            }
            return sum;
        }

        //Dixon-Price  50D  -10,10   min = 0
        private double Function13(double[] XD)
        {
            double sum = 0;
            for (int i = 1; i < XD.Length; i++)
            {
                double p1 = 2 * Math.Pow(XD[i], 2);
                sum += i * Math.Pow(p1 - XD[i - 1], 2);
            }
            sum += Math.Pow(XD[0] - 1, 2);
            return sum;
        }

        //Powell Singular  D50 range -4,5 min=0
        private double Function14(double[] XD)
        {
            double sum = 0;

            for (int i = 0; i < XD.Length / 4; i++)
            {
                double p1 = (XD[4 * (i + 1) - 3 - 1] + 10 * XD[4 * (i + 1) - 2 - 1]);
                double p2 = 5 * (XD[4 * (i + 1) - 1 - 1] - XD[4 * (i + 1) - 1]);
                double p3 = XD[4 * (i + 1) - 2 - 1] - XD[4 * (i + 1) - 1 - 1];
                double p4 = 10 * (XD[4 * (i + 1) - 3 - 1] - XD[4 * (i + 1) - 1]);
                sum = sum + (Math.Pow(p1, 2) + Math.Pow(p2, 2) + Math.Pow(p3, 4) + Math.Pow(p4, 4));
            }
            return sum;
        }
        //Zakharov D50  range -5,10 min = 0
        private double Function15(double[] XD)
        {
            double sum = 0;
            double p1;
            double p2;
            p1 = p2 = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(XD[i], 2);
                p2 += (0.5 * (i + 1) * XD[i]);
            }

            sum = p1 + Math.Pow(p2, 2) + Math.Pow(p2, 4);
            return sum;
        }

        // Xin-She Yang 3 50D -20,20 min = -1
        private double Function16(double[] XD)
        {
            double sum = 0;
            double beta = 15;
            double m = 5;
            double p1 = 0;
            double p2 = 0;
            double p3 = 1;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(XD[i] / beta, 2 * m);
                p2 += Math.Pow(XD[i], 2);
                p3 *= Math.Pow(Math.Cos(XD[i]), 2);
            }

            sum = (Math.Exp(-p1) - 2 * Math.Exp(-p2)) * p3;
            return sum;
        }

        //Perm 0, D, Beta 5D, range -D, D,  min = 0, X= {1, 1/2, ..., 1/n}
        private double Function17(double[] XD)
        {
            double sum = 0;
            double p1 = 0;
            double beta = 10;
            for (int i = 0; i < XD.Length; i++)
            {
                for (int j = 0; j < XD.Length; j++)
                {
                    //p1 += (j + 1 + beta) * (Math.Pow(XD[j], i + 1) - (1 / Math.Pow(j + 1, i + 1)));
                    double jB = j + 1 + beta;
                    double xp = Math.Pow(XD[j], i + 1);
                    double jP = Math.Pow(j + 1, i + 1);
                    p1 += jB * (xp - (1 / jP));
                }
                sum += Math.Pow(p1, 2);
            }
            return sum;
        }
        // Three-Hump Camel 2D, -5,5, min = 0
        private double Function18(double[] XD)
        {
            double sum = 0.0f;
            sum = 2 * Math.Pow(XD[0], 2) - 1.05 * Math.Pow(XD[0], 4) +
                   Math.Pow(XD[0], 6) / 6 +
                   XD[0] * XD[1] + Math.Pow(XD[1], 2);
            return sum;
        }
        //Beale 2D , range -4.5, 4.5 min = 0
        private double Function19(double[] XD)
        {
            double sum = 0.0f;
            double p2 = Math.Pow(XD[1], 2);
            double p3 = Math.Pow(XD[1], 3);
            sum = Math.Pow((1.5 - XD[0] + (XD[0] * XD[1])), 2) +
                  Math.Pow((2.25 - XD[0] + (XD[0] * p2)), 2) +
                  Math.Pow((2.625 - XD[0] + (XD[0] * p3)), 2);

            return sum;
        }
        //Booth  2D -10,10 min = 0
        private double Function20(double[] XD)
        {
            double sum = 0;
            sum = Math.Pow((XD[0] + (2 * XD[1]) - 7), 2) + Math.Pow(((2 * XD[0]) + XD[1] - 5), 2);
            return sum;
        }

        //Brent  2D -10,10   min = 0
        private double Function21(double[] XD)
        {
            double sum = 0;
            double p1 = Math.Pow(XD[0], 2);
            double p2 = Math.Pow(XD[1], 2);
            sum = Math.Pow(XD[0] + 10, 2) + Math.Pow(XD[1] + 10, 2) + Math.Exp(-p1 - p2);
            return sum;
        }
        //Matyas, 2D, range -10,10 min = 0
        private double Function22(double[] XD)
        {
            double a = Math.Pow(XD[0], 2) + Math.Pow(XD[1], 2);

            double sum = (0.26 * a) - (0.48 * XD[0] * XD[1]);
            return sum;

        }
        //Schaffer N.4  2D,  -100,100   min = 0.29257
        private double Function23(double[] XD)
        {

            double sum = 0;
            double sin = Math.Sin(Math.Abs(XD[0] * XD[0] - XD[1] * XD[1]));
            double cos = Math.Cos(sin);
            double pow = Math.Pow(cos, 2);
            double up = pow - 0.5;

            double p = 1 + 0.001 * (XD[0] * XD[0] + XD[1] * XD[1]);
            double dawn = Math.Pow(p, 2);

            sum = 0.5 + (up / dawn);
            return sum;
        }
        //Wayburn seader 3,  20 D, range [-500, 500],  min = 9.10
        private double Function24(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = Math.Pow(XD[0], 2);
            double p2 = Math.Pow(XD[0], 3);
            double p3 = Math.Pow(XD[0] - 4, 2) + Math.Pow(XD[1] - 5, 2) - 4;
            sum = (2 * p2 / 3) - (8 * p1) + (33 * XD[0]) - (XD[0] * XD[1]) + 5 + Math.Pow(p3, 2);

            return sum;
        }
        //Leon,  2D,  -1.2, 1.2,    min = 0
        private double Function25(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = Math.Pow((XD[1] - Math.Pow(XD[0], 3)), 2);
            double p2 = Math.Pow(1 - XD[0], 2);
            sum = 100 * p1 + p2;
            return sum;
        }


        /**********************************
        *
        * Multimodal Bench mark functions 
        *
        **********************************/
        //Schwefel 2.26 50 D,  -500, 500 min = 0
        private double Function26(double[] XD)
        {
            double sum = 0;
            double z = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                double abs = Math.Abs(XD[i]);
                double sqrt = Math.Sqrt(abs);
                double sin = Math.Sin(sqrt);
                z += ((XD[i] * sin));
            }
            sum = 418.9829 * XD.Length - z;
            return sum;
        }

        //Rastrigin  50 D   -5.12, 5.12  min = 0
        private double Function27(double[] XD)
        {
            double sum = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                sum += ((XD[i] * XD[i]) - 10 * Math.Cos(2 * XD[i] * Math.PI) + 10);
            }
            return sum;
        }
        //Periodic 50 D   -10,10   min = 0.9
        private double Function28(double[] XD)
        {
            double sum = 0;
            double p1 = 0;
            double p2 = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(Math.Sin(XD[i]), 2);
                p2 += Math.Pow(XD[i], 2);
            }
            sum = 1 + p1 - (0.1 * Math.Exp(-p2));
            return sum;
        }

        // Qing 50 D,   -500,500   min = 0   X= -+ Squrt(i)
        private double Function29(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                double p = Math.Pow(XD[i], 2);
                sum += Math.Pow(p - (i + 1), 2);
            }
            return sum;
        }
        //Alpine N.1 50 D,  -10, 10   min = 0
        private double Function30(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
                sum += Math.Abs(XD[i] * Math.Sin(XD[i]) + 0.1 * XD[i]);
            return sum;
        }

        //Xin-She Yang N.5  50 D, -5, 5 min = 0  X= 0, 0...0
        private double Function31(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                double abs = Math.Abs(XD[i]);
                double pow = Math.Pow(abs, i + 1);
                sum += generater.NextDouble() * pow;
            }

            return sum;
        }

        //Ackley N.1 50 D, -32, 32   min 0   X = 0, 0 ...0
        private double Function32(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(XD[i], 2);
                p2 += Math.Cos(2 * Math.PI * XD[i]);
            }
            double sqrt = Math.Sqrt(p1 / XD.Length);
            sum = -20 * Math.Exp(-0.2 * sqrt) - Math.Exp(p2 / XD.Length) + 20 + Math.Exp(1);
            return sum;
        }

        // Trignometric 2   50 D,  -500, 500 min = 1 X= 
        private double Function33(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                double pow1 = 7 * Math.Pow(XD[i] - 0.9, 2);
                double p1 = 8 * Math.Pow(Math.Sin(pow1), 2);
                double pow2 = 14 * Math.Pow(XD[i] - 0.9, 2);
                double p2 = 6 * Math.Pow(Math.Sin(pow2), 2);
                double p3 = Math.Pow(XD[i] - 0.9, 2);
                sum += (p1 + p2 + p3);
            }
            return sum;
        }

        //Salomon  50 D,  -100, 100  min = 0   X=0, 0,   0
        private double Function34(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p = 0;
            for (int i = 0; i < XD.Length; i++)
                p += Math.Pow(XD[i], 2);
            double sqrt = Math.Sqrt(p);
            sum = 1 - Math.Cos(2 * Math.PI * sqrt) + (0.1 * sqrt);
            return sum;
        }
        // Styblinski-Tang 50 D, -5, 5   min = -39.16599*n
        private double Function35(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            for (int i = 0; i < XD.Length; i++)
            {
                double p1 = Math.Pow(XD[i], 4);
                double p2 = 16 * Math.Pow(XD[i], 2);
                double p3 = 5 * XD[i];
                sum += p1 - p2 + p3;
            }
            sum = sum * 0.5;
            return sum;
        }
        // Griewank  50 D,  -100, 100,  min = 0  X = 0, 0...0
        private double Function36(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 1;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(XD[i], 2);
                double sqrt = Math.Sqrt(i + 1);
                double cos = Math.Cos(XD[i] / sqrt);
                p2 *= cos;
            }
            sum = (p1 / 4000) - p2 + 1;
            return sum;
        }
        //Xin-She Yang N.4  50 D, -10, 10   min = -1
        private double Function37(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 0;
            double p3 = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Pow(Math.Sin(XD[i]), 2);
                p2 += Math.Pow(XD[i], 2);
                double sqrt = Math.Sqrt(Math.Abs(XD[i]));
                double sin = Math.Pow(Math.Sin(sqrt), 2);
                p3 += sin;
            }
            sum = (p1 - Math.Exp(-p2)) * Math.Exp(-p3);
            return sum;
        }
        //Xin-She Yang N.3 50 D   -2Pi, 2Pi  min = 0, X = 0, 0,..0
        private double Function38(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 0;
            for (int i = 0; i < XD.Length; i++)
            {
                p1 += Math.Abs(XD[i]);
                p2 += Math.Sin(Math.Pow(XD[i], 2));
            }
            sum = p1 * Math.Exp(-p2);
            return sum;
        }

        // Generalized Penalized  50 D, -50, 50  min = 0, X = 1, 1,.. 1
        private double Function39(double[] XD)
        {
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 0;
            double p3 = 0;
            double p4 = 0;
            double sum1 = 0;
            double sum2 = 0;
            p3 = Math.Pow(XD[XD.Length - 1] - 1, 2);
            p4 = 1 + Math.Pow(Math.Sin(2 * Math.PI * XD[XD.Length - 1]), 2);

            for (int i = 0; i < XD.Length - 1; i++)
            {
                p1 = Math.Pow(XD[i] - 1, 2);
                p2 = 1 + Math.Pow(Math.Sin(3 * Math.PI * XD[i + 1]), 2);
                //p3 = Math.Pow(XD[XD.Length - 1] - 1, 2);
                //p4 = 1+Math.Pow(Math.Sin(2*Math.PI*XD[XD.Length-1]),2);
                sum1 += p1 * p2;
            }
            sum1 += p3 * p4;

            for (int i = 0; i < XD.Length; i++)
            {
                double r = 0;
                if (XD[i] > 5)
                    r = 100 * Math.Pow(XD[i] - 5, 4);
                else if (XD[i] < -5)
                    r = 100 * Math.Pow(-XD[i] - 5, 4);
                else
                    r = 0;
                sum2 += r;
            }
            sum = 0.1 * (Math.Pow(Math.Sin(3 * Math.PI * XD[0]), 2) + sum1) + sum2;
            return sum;
        }
        //Generalized Penalized N.1  50 D, -50, 50  min = 0,  X = -1
        private double Function40(double[] XD)
        {
            double sum = 0.0f;
            double p1 = 0;
            double p2 = 0;
            double sum1 = 0;
            double sum2 = 0;
            double y1 = 0;
            double y2 = 0;
            for (int i = 0; i < XD.Length - 1; i++)
            {
                y1 = 1 + (XD[i] + 1) / 4;
                y2 = 1 + (XD[i + 1] + 1) / 4;
                p1 = Math.Pow(y1 - 1, 2);
                p2 = 1 + 10 * Math.Pow(Math.Sin(Math.PI * y2), 2);
                sum1 += p1 * p2;
            }

            for (int i = 0; i < XD.Length; i++)
            {
                double r = 0;
                if (XD[i] > 10)
                    r = 100 * Math.Pow(XD[i] - 10, 4);
                else if (XD[i] < -10)
                    r = 100 * Math.Pow(-XD[i] - 10, 4);
                else
                    r = 0;
                sum2 += r;
            }
            double yy = 1 + (XD[0] + 1) / 4;
            double yn = 1 + (XD[XD.Length - 1] + 1) / 4;
            sum = (Math.PI / XD.Length) * ((10 * Math.Pow(Math.Sin(Math.PI * yy), 2) + sum1 + Math.Pow(yn - 1, 2) + sum2));
            return sum;
        }
        //Egg crate 2D, -5,5  min = 0 
        private double Function41(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            sum = Math.Pow(XD[0], 2) + Math.Pow(XD[01], 2) +
                  25 * (Math.Pow(Math.Sin(XD[0]), 2) + Math.Pow(Math.Sin(XD[1]), 2));
            return sum;
        }

        //Ackley N.3  2D, -32,32  min = -195.629028238419  X= +- 0.682584587366, -0.3607532551372 
        private double Function42(double[] XD)
        {
            //Ok

            double sum = 0.0f;
            double p1 = -0.2 * Math.Sqrt((XD[0] * XD[0]) + (XD[1] * XD[1]));
            double p2 = Math.Cos(3 * XD[0]) + Math.Sin(3 * XD[1]);
            sum = (-200 * Math.Exp(p1)) + (5 * Math.Exp(p2));
            return sum;
            //return -200*np.exp(-0.2*np.sqrt(x**2 + y**2)) + 
            //5*np.exp(np.cos(3*x)+np.sin(3*y))
        }
        // Adjiman 2D  -1, 2  min = -2.02181 X= 0,0 
        private double Function43(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            sum = Math.Cos(XD[0]) * Math.Sin(XD[1]) - (XD[0] / (Math.Pow(XD[1], 2) + 1));
            return sum;
        }
        //Bird 2D   -2PI,2PI    min = -106.7645  X= -1.58214, -3.13024
        private double Function44(double[] XD)
        {
            //Ok
            double sum = 0.0f;
            double p1 = Math.Pow(1 - Math.Cos(XD[1]), 2);
            double p2 = Math.Pow(1 - Math.Sin(XD[0]), 2);
            double p3 = Math.Pow(XD[0] - XD[1], 2);
            sum = Math.Sin(XD[0]) * Math.Exp(p1) + Math.Cos(XD[1]) * Math.Exp(p2) + p3;
            return sum;
        }

        // Camel Six Hump,   2D,  -5, 5   min = -1.0316  (x1,x2)=(-0.0898,0.7126), (0.0898,-0.7126)
        private double Function45(double[] XD)
        {

            double sum = 0;
            sum = 4 * Math.Pow(XD[0], 2) - 2.1 * Math.Pow(XD[0], 4) + Math.Pow(XD[0], 6) / 3 + XD[0] * XD[1] - 4 * Math.Pow(XD[1], 2) + 4 * Math.Pow(XD[1], 4);
            return sum;
        }

        // Branin RCOS  2D ,  x1[-5, 10], x2[0,15]  min = 0.397887   X= +- (pi, 12.275)
        private double Function46(double[] XD)
        {

            double sum = 0;
            double p1 = XD[1] - (5.1 * Math.Pow(XD[0], 2) / (4 * Math.PI * Math.PI)) + (5 * XD[0] / Math.PI) - 6;
            p1 = Math.Pow(p1, 2);
            double p2 = 10 * (1 - 1 / (8 * Math.PI)) * Math.Cos(XD[0]);
            sum = p1 + p2 + 10;
            return sum;
        }
        //Hartman 3  3D,  0, 1  min = -3.862782  X = 0.114614, 0.555649, 0.852547
        private double Function47(double[] XD)
        {
            double[] a = { 1, 1.2, 3, 3.2 };
            double[,] A = { {3, 10, 30},
                            {0.1, 10, 35},
                            {3, 10, 30},
                            {0.1, 10, 35}
                        };

            double[,] p = {{ 3689, 1170, 2673},
                          { 4699, 4387, 7470},
                          { 1091, 8732, 5547},
                          {381, 5743, 8828} };
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 3; j++)
                    p[i, j] *= Math.Pow(10, -4);
            double sum = 0;
            double sum1 = 0;
            for (int i = 0; i < 4; i++)
            {
                sum1 = 0;
                for (int j = 0; j < 3; j++)
                    sum1 += A[i, j] * Math.Pow(XD[j] - p[i, j], 2);
                sum += a[i] * Math.Exp(-sum1);
            }
            sum *= -1;
            return sum;
        }
        //Hartman 6  6D,  0, 1  min = -3.3237  X = 0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573
        private double Function48(double[] XD)
        {
            double[] a = { 1, 1.2, 3, 3.2 };
            double[,] A = { {10, 3, 17, 3.5, 1.7, 8},
                            {0.05, 10, 17, 0.1, 8, 14},
                            {3, 3.5, 1.7, 10, 17, 8},
                            {17, 8, 0.05, 10, 0.1, 14}
                        };

            double[,] p = {{ 1312, 1696, 5569, 124, 8283, 5886},
                          { 2329, 4135, 8307, 3736, 1004, 9991},
                          { 2348, 1451, 3522, 2883, 3047, 6650},
                          {4047, 8828, 8732, 5743, 1091, 381} };
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 6; j++)
                    p[i, j] *= Math.Pow(10, -4);
            double sum = 0;
            double sum1 = 0;
            for (int i = 0; i < 4; i++)
            {
                sum1 = 0;
                for (int j = 0; j < 6; j++)
                    sum1 += (A[i, j] * Math.Pow(XD[j] - p[i, j], 2));
                sum += a[i] * Math.Exp(-sum1);
            }
            sum *= -1;
            return sum;
        }
        // Cross-in-tray  2D  -10,10 min = -2.06261218  X = +-1.34940668535334, +-1.349406608602064
        private double Function49(double[] XD)
        {

            double sum = 0;
            double p1 = Math.Abs(100 - (Math.Sqrt(XD[0] * XD[0] + XD[1] * XD[1])) / Math.PI);
            double exp = Math.Exp(p1);
            double sinx = Math.Sin(XD[0]);
            double siny = Math.Sin(XD[1]);
            double abs = Math.Abs(sinx * siny * exp) + 1;
            double pow = Math.Pow(abs, 0.1);
            sum = -0.0001 * pow;
            return sum;
        }

        //Bartels Conn 2D,  -500, 500  min = 1,  X = 0, 0
        private double Function50(double[] XD)
        {

            double sum = 0;
            double p1 = Math.Abs(XD[0] * XD[0] + XD[1] * XD[1] + XD[0] * XD[0]);
            double p2 = Math.Abs(Math.Sin(XD[0]));
            double p3 = Math.Abs(Math.Cos(XD[1])); ;

            sum = p1 + p2 + p3;
            return sum;
        }
        private double GN(Random random, double mean, double stddev)
        {
            // The method requires sampling from a uniform random of (0,1]
            // but Random.NextDouble() returns a sample of [0,1).
            double x1 = 1 - random.NextDouble();
            double x2 = 1 - random.NextDouble();

            double y1 = Math.Sqrt(-2.0 * Math.Log(x1)) * Math.Cos(2.0 * Math.PI * x2);
            return (y1 * stddev + mean);
            //return y1 ;
        }

    }//end class
    internal class Program
    {
        static void Main(string[] args)
        {
            Console.Clear();
            Console.WriteLine("Please wait while running....");

            int Pubulaionsize = 30;
            int Iteration = 500;
            int Runs = 10;
            /****************************************************
             * This code is using (Sphere function) only
             * Sphere 50 D,  range = [-100, 100] 
             *   Min = 0,     x = 0,0....0
             *    
             * ***************************************************/
            double upper = 100;// upper of Function 1 (Sphere function)
            double lower = -100;// Lower of Function 1
            double global = 0;// Golobal optimum of Function 1
            int dimention = 50;// 50 dimention
            double stopCondition = 0;

            double[] BP = null;
            List<double> BestSoFar = new List<double>();
            //List<double[]> LastGenerations = new List<double[]>();
            List<double> NOFE = new List<double>();
            int succes = 0;

            CHE ct = null;
            for (int i = 0; i < Runs; i++)
            {

                ct = new CHE(dimention, Pubulaionsize, Iteration, upper, lower, global, stopCondition);
                ct.run();
                BestSoFar.Add(ct.BestFintness());
                NOFE.Add(ct.GetNOFE());
                if (ct.Succes)
                {
                    succes++;

                }
                string sbestFitness = "Run: " + i.ToString() + "  Best fitness = " + ct.BestFintness().ToString();
                Console.WriteLine(sbestFitness);
                Console.WriteLine("NOFA: " + ct.GetNOFE().ToString());

                BP = new double[ct.GeBestPoints().Length];
                BP = ct.GeBestPoints();
                string s = "Pest points: {";
                for (int j = 0; j < BP.Length; j++)
                    s += BP[j].ToString() + " , ";

                s += "}";
                Console.WriteLine(s);
                Console.WriteLine();

            }

            Console.WriteLine("-------------------------------------------");
            Console.WriteLine();

            double AvBestSoFar = 0;
            double best_so_far = BestSoFar.Min();
            if (BestSoFar.Count > -1)
                AvBestSoFar = BestSoFar.Average();

            BestSoFar.Sort();
            double Median = BestSoFar[BestSoFar.Count / 2];
            double sumNOFE = 0;
            double meanNOFE = 0;
            double sdNOFE = 0;
            if (NOFE.Count >= 0)
            {
                for (int i = 0; i < NOFE.Count; i++)
                    sumNOFE += NOFE[i];
                meanNOFE = sumNOFE / NOFE.Count;

                double x = 0;
                for (int i = 0; i < NOFE.Count; i++)
                    x = x + Math.Pow((NOFE[i] - meanNOFE), 2);
                double sdSum = x / (double)(NOFE.Count);
                sdNOFE = Math.Sqrt(sdSum);
            }

            Console.WriteLine("-------------------------------");
            Console.WriteLine();
            float suc = (float)succes / Runs;
            Console.WriteLine("Succes Rate = " + suc.ToString());
            Console.WriteLine(" best so far fitness = " + best_so_far.ToString());
            Console.WriteLine("Average best so far fitness = " + AvBestSoFar.ToString());
            Console.WriteLine("Mediam best so far fitness = " + Median.ToString());
            Console.WriteLine("Mean NOFE = " + meanNOFE.ToString());
            Console.WriteLine("NOFE Standard Divetion = " + sdNOFE.ToString());

            Console.ReadLine();


        }

    }
}