#include <iostream>
#include <fstream>
#include <vector>
#include "Timer.h"

#define DF2

#ifdef DF1
typedef float type, sum;
#endif
#ifdef DF2
typedef double type, sum;
#endif
#ifdef DF3
typedef float type;
typedef double sum;
#endif

struct matrix
{
public:
   matrix()
   {
      this->N_ = 0;
      this->maxiter_ = 0;
      this->epsilon_ = 0;
      this->ig_ = std::vector<int>();
      this->jg_ = std::vector<int>();
      this->ggl_ = std::vector<type>();
      this->ggu_ = std::vector<type>();
      this->di_ = std::vector<type>();
      this->pr_ = std::vector<type>();
      this->x_ = std::vector<type>();
   }
   matrix(int N,
      int maxiter,
      type epsilon,
      std::vector<int>  ig,
      std::vector<int> jg,
      std::vector<type> ggl,
      std::vector<type> ggu,
      std::vector<type> di,
      std::vector<type> pr,
      std::vector<type> x
   )
   {
      N_ = N;
      maxiter_ = maxiter;
      epsilon_ = epsilon;
      ig_ = ig;
      jg_ = jg;
      ggl_ = ggl;
      ggu_ = ggu;
      di_ = di;
      pr_ = pr;
      x_ = x;
   }

   void load_data()
   {
      std::ifstream di,ggl,ggu,ig,jg,kuslau,pr,x, xr;
      
      kuslau.open("kuslau.txt");
      kuslau >> N_ >> maxiter_ >> epsilon_;
      kuslau.close();

      ig.open("ig.txt");
      ig_.resize(N_ + 1);
      for (int i = 0; i < ig_.size(); ++i)
      {
         ig >> ig_[i];
         ig_[i]--;
      }
      ig.close();

      jg.open("jg.txt");
      jg_.resize(ig_.back());
      for (int i = 0; i < jg_.size(); ++i)
      {
         jg >> jg_[i];
         jg_[i]--;
      }
      jg.close();

      ggl.open("ggl.txt");
      ggl_.resize(ig_.back());
      for (int i = 0; i < ggl_.size(); ++i)
      {
         ggl >> ggl_[i];
      }
      ggl.close();

      ggu.open("ggu.txt");
      ggu_.resize(ig_.back());
      for (int i = 0; i < ggu_.size(); ++i)
      {
         ggu >> ggu_[i];
      }
      ggu.close();

      di.open("di.txt");
      di_.resize(N_);
      for (int i = 0; i < N_; ++i)
      {
         di >> di_[i];
      }
      di.close();

      pr.open("pr.txt");
      pr_.resize(N_);
      for (int i = 0; i < N_; ++i)
      {
         pr >> pr_[i];
      }
      pr.close();

      x.open("x.txt");
      x_.resize(N_);
      for (int i = 0; i < N_; ++i)
      {
         x >> x_[i];
      }
      x.close();

      //print_(di_, ggu_, ggl_);
   }

   void print_(std::vector<type>& di, std::vector<type>& ggu, std::vector<type>& ggl)
   {
      std::vector <std::vector<type>> mat;
      mat.resize(N_);
      for (int i = 0; i < mat.size(); i++)
      {
         mat[i].resize(N_);
      }

      for (int i = 0; i < mat.size(); i++)
      {
         mat[i][i] = di[i];
         for (int j = ig_[i]; j < ig_[i + 1]; j++)
         {
            mat[i][jg_[j]] = ggl[j];
            mat[jg_[j]][i] = ggu[j];
         }
      }

      for (int i = 0; i < mat.size(); i++)
      {
         for (int j = 0; j < mat[i].size(); j++)
            std::cout << mat[i][j] << " ";
         std::cout << std::endl;
      }
   }

   void LUsq()
   {
      std::vector<type> di_f;
      std::vector<type> ggl_f;
      std::vector<type> ggu_f;

      di_f = di_;
      ggl_f = ggl_;
      ggu_f = ggu_;
      
      for (int i = 0; i < N_; i++)
      {
         type sumdi = 0.0;		

         int i0 = ig_[i];		
         int i1 = ig_[i + 1];	

         
         for (int k = i0; k < i1; k++)
         {
            int j = jg_[k];		
            int j0 = ig_[j];		
                                 
            int j1 = ig_[j + 1];	
                                 

            int ik = i0;			
            int kj = j0;			

            type suml = 0.0;		
            type sumu = 0.0;		

            while (ik < k)
            {
               
               if (jg_[ik] == jg_[kj])
               {
                 
                  suml += ggl_f[ik] * ggu_f[kj];
                  sumu += ggu_f[ik] * ggl_f[kj];
                  ik++;
                  kj++;
               }
              
               else
                  jg_[ik] > jg_[kj] ? kj++ : ik++;
            }

            ggl_f[k] = (ggl_f[k] - suml) / di_f[j];
            ggu_f[k] = (ggu_f[k] - sumu) / di_f[j];
            sumdi += ggl_f[k] * ggu_f[k];
         }

         di_f[i] = sqrt(di_f[i]-sumdi);
      }
 //     print_(di_f, ggu_f, ggl_f);
      LoS_precond(di_f, ggu_f, ggl_f);
   }

   std::vector<type> Mult(std::vector<type>& v)
   {
      std::vector<type> res(v.size());
      for (int i = 0; i < v.size(); i++)
      {
         res[i] = di_[i] * v[i];
         for (int j = ig_[i]; j < ig_[i + 1]; j++)
         {
            res[i] += ggl_[j] * v[jg_[j]];
            res[jg_[j]] += ggu_[j] * v[i];
         }
      }
      return res;
   }

   void LoS()
   {
      type error = 100;
      int k = 1;
      std::vector<type> r;
      std::vector<type> buf = Mult(x_);
      r.resize(N_);
      for (int i = 0; i < N_; i++)
      {
         r[i] = pr_[i] - buf[i];
      }
      error = scalar_prod(r, r);
      std::vector<type> z = r;
      std::vector<type> p = Mult(z);
      while (error>epsilon_ && k < maxiter_)
      {
         type pp = scalar_prod(p, p);
         type alpha = scalar_prod(p, r) / pp;
         error -= alpha * alpha * pp;
         for (int i = 0; i < N_; i++)
         {
            x_[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         std::vector<type> buf = Mult(r);
         type betta = -(scalar_prod(p, buf) / pp);
         for (int i = 0; i < N_; i++)
         {
            z[i] = r[i] + betta * z[i];
            p[i] = buf[i] + betta * p[i];
         }
         //print_console(k, error);
         k++;
      }
      std::cout << "k:" << k<< "\nerror: " << error << "\n";
      print_file(k, error);
   }

   void LoS_precond(std::vector<type>& di_f, std::vector<type>& ggu_f, std::vector<type>& ggl_f)
   {
      type error = 100;
      int k = 1;
      std::vector<type> buf = Mult(x_);
      for (int i = 0; i < N_; i++)
      {
         buf[i] = pr_[i] - buf[i];
      }
      std::vector<type> r = LUDirect(buf, di_f, ggl_f);
      error = scalar_prod(r, r);
      std::vector<type> z = LUReverse(r, di_f, ggu_f);
      buf = Mult(z);
      std::vector<type> p = LUDirect(buf, di_f, ggl_f);;
      while (error > epsilon_ && k < maxiter_)
      {
         type pp = scalar_prod(p, p);
         type pr = scalar_prod(p, r);
         type alpha = pr / pp;
         error -= alpha * alpha * pp;
         for (int i = 0; i < N_; i++)
         {
            x_[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         std::vector<type> Ur = LUReverse(r, di_f, ggu_f);
         buf = Mult(Ur);
         buf = LUDirect(buf, di_f, ggl_f);
         type betta = -(scalar_prod(p, buf) / pp);
         for (int i = 0; i < N_; i++)
         {
            z[i] = Ur[i] + betta * z[i];
            p[i] = buf[i] + betta * p[i];
         }
         //print_console(k, error);
         k++;
      print_file(k, error);
      }
      std::cout << "k:" << k << "\nerror: " << error << "\n";
   };

   std::vector<type> LUDirect(const std::vector<type>& b, std::vector<type>& di_f, std::vector<type>& ggl_f)
   {
      std::vector<type> res = b;

      for (size_t i = 0; i < res.size(); i++)
      {
         type sum = 0.0;
         for (size_t j = ig_[i]; j < ig_[i + 1]; j++)
            sum += ggl_f[j] * res[jg_[j]];
         res[i] -= sum;
         res[i] /= di_f[i];
      }
      return res;
   }

   std::vector<type> LUReverse(const std::vector<type>& b, std::vector<type>& di_f, std::vector<type>& ggu_f)
   {
      std::vector<type> res = b;

      for (int i = N_ - 1; i >= 0; i--)
      {
         res[i] /= di_f[i];
         for (size_t j = ig_[i]; j < ig_[i + 1]; j++)
            res[jg_[j]] -= ggu_f[j] * res[i];
      }
      return res;
   }

   type scalar_prod(std::vector<type>& x, std::vector<type>& y)
   {
      type res = 0.0;
      if (x.size() == y.size())
      {
         for (int i = 0; i < y.size(); i++)
         {
            res += x[i] * y[i];
         }
         return res;
      }
      else
      {
         std::cout << "Error!";
         return res;
      }
   }

   void print_console(int k, type error)
   {
      std::cout << "k: " << k << std::endl;
      std::cout << "nevizka: " << error << std::endl;
      for (int i = 0; i < N_; i++)
      {
         std::cout << x_[i] << std::endl;
      }
   }

   void print_file(int k, type error)
   {
      std::fstream out;
      out.open("output.txt", std::ios_base::app);
      out << "k: " << k << std::endl;
      out << "nevizka: " << error << std::endl;
      for (int i = 0; i < N_; i++)
      {
         out << x_[i] << std::endl;
      }
   }
   
   std::vector<type> diagFact(std::vector<type> di)
   {
      std::vector<type> di_f = di;
      for (size_t i = 0; i < N_; i++)
         di_f[i] = 1.0 / sqrt(di_f[i]);
      return di_f;
   }

   std::vector<type> multD(const std::vector<type> v, std::vector<type> di_f)
   {
      std::vector<type> res(N_);
      for (size_t i = 0; i < N_; i++)
         res[i] = di_f[i] * v[i];
      return res;
   }

   void LoS_diag()
   {
      std::vector<type> di_f = diagFact(di_);
      type error = 100;
      int k = 1;
      std::vector<type> buf = Mult(x_);
      for (int i = 0; i < N_; i++)
      {
         buf[i] = pr_[i] - buf[i];
      }
      std::vector<type> r = multD(buf, di_f);
      error = scalar_prod(r, r);
      std::vector<type> z = multD(r,di_f);
      buf = Mult(z);
      std::vector<type> p = multD(buf, di_f);;
      while (error > epsilon_ && k < maxiter_)
      {
         type pp = scalar_prod(p, p);
         type pr = scalar_prod(p, r);
         type alpha = pr / pp;
         error -= alpha * alpha * pp;
         for (int i = 0; i < N_; i++)
         {
            x_[i] += alpha * z[i];
            r[i] -= alpha * p[i];
         }
         std::vector<type> Ur = multD(r, di_f);
         buf = Mult(Ur);
         buf = multD(buf, di_f);
         type betta = -(scalar_prod(p, buf) / pp);
         for (int i = 0; i < N_; i++)
         {
            z[i] = Ur[i] + betta * z[i];
            p[i] = buf[i] + betta * p[i];
         }
         //print_console(k, error);
         k++;
      }
      std::cout << "k:" << k << "\nerror: " << error << "\n";
      print_file(k, error);
   };

private:
   int N_;
   int maxiter_;
   type epsilon_;
   std::vector<int> ig_;
   std::vector<int> jg_;
   std::vector<type> ggl_;
   std::vector<type> ggu_;
   std::vector<type> di_;
   std::vector<type> pr_;
   std::vector<type> x_;
};

int main()
{
   matrix Mx;
   Mx.load_data();
   Timer timer;
   Mx.LUsq();
   //Mx.LoS_diag();
  //Mx.LoS();
   return 0;
}

