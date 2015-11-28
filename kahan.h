// This is taken from the OpenSAGE project archived on Google Code at
// http://opensage.googlecode.com/svn/core_code/trunk/src/c++/include/numerics/kahan.h

#ifndef KAHAN_H
#define KAHAN_H

// Kahan Summation (also known as compensated summation) template
//
// This technique is due to Dr. William Kahan, University of California,
//   Berkeley
//
// Ref: 
//   "What Every Computer Scientist Should Know About Floating Point
//    Arithmetic," by David Goldberg, published in the March, 1991 issue of
//    Computing Surveys. Copyright 1991, ACM
//
// Theorem 8 (Kahan Summation Formula) on page 203 shows:
//  error in iterative summation can be reduced to 2 * epsilon,
//  whereas naive summation can produce n * epsilon.
//

// NOTE:  This is the most general template that assumes a very large or
//        dynamically allocated class T.  So large that temporaries are
//        pre-allocated in the class. See the partial specializations below
//        for small floating point values (double, float).

template <class T> 
class KahanAdder
{
public:
      KahanAdder()                      : my_sum(0), my_carry(0)  {}
      KahanAdder(const T& e)            : my_sum(e), my_carry(0)  {}
      KahanAdder(const T& e, const T& c): my_sum(e), my_carry(c)  {}

      // return the working sum.
      const T &sum() const { return my_sum; }
      // return the current carry.
      const T &carry() const { return my_carry; }

      operator const T&() const { return my_sum; }

      template<class S>
      KahanAdder &operator =  (const S &g) 
      { my_sum = g; my_carry = 0; return *this; }

      template<class S>
      KahanAdder &operator += (const S &g)  { add(g); return *this; }
      template<class S>
      KahanAdder &operator -= (const S &g)  { add(-g); return *this; }
      template<class S>
      KahanAdder &operator *= (const S &g)  { *this = ((T)*this) * g; return *this; }
      template<class S>
      KahanAdder &operator /= (const S &g)  { *this = ((T)*this) / g; return *this; }

private:
      T my_sum;   // The current working sum.
      T my_carry; // The carry from the previous operation

      // Temporaries are members due to construction overhead
      T my_temp;
      T my_y;

      // After each add operation, the carry will contain the additional
      // bits that could be left out from the previous operation.
      void add(const T &datum) 
      {
          my_y     = datum - my_carry;
          my_temp  = my_sum + my_y;
          my_carry = (my_temp - my_sum) - my_y;
          my_sum   = my_temp;
      }
};

///////////////////////////////////////////////////////////////////////
//                 SPECILIZATIONS OF KahanAdder<T>                   //
//-------------------------------------------------------------------//
//        Lets be more efficient for T = {double, float} by          //
//  removing the temporaries that are stored in the class.  These    //
//     are only a win when the floating point type is huge and/or    //
//       requires dynamic allocation to create a new instance.       //
///////////////////////////////////////////////////////////////////////

template <> 
class KahanAdder<double>
{
public:
      KahanAdder()                  : my_sum(0), my_carry(0)  {}
      KahanAdder(double e)          : my_sum(e), my_carry(0)  {}
      KahanAdder(double e, double c): my_sum(e), my_carry(c)  {}

      // return the working sum.
      const double &sum() const { return my_sum; }
      // return the current carry.
      const double &carry() const { return my_carry; }

      operator const double&() const { return my_sum; }

      template<class S>
      KahanAdder<double> &operator =  (const S &g) 
      { my_sum = g; my_carry = 0; return *this; }

      template<class S>
      KahanAdder<double> &operator += (const S &g)  { add(g); return *this; }
      template<class S>
      KahanAdder<double> &operator -= (const S &g)  { add(-g); return *this; }
      template<class S>
      KahanAdder<double> &operator *= (const S &g)  { *this = ((double)*this) * g; return *this; }
      template<class S>
      KahanAdder<double> &operator /= (const S &g)  { *this = ((double)*this) / g; return *this; }

private:
      double my_sum;   // The current working sum.
      double my_carry; // The carry from the previous operation

      // After each add operation, the carry will contain the additional
      // bits that could be left out from the previous operation.
      void add(const double &datum) 
      {
          double my_y     = datum - my_carry;
          double my_temp  = my_sum + my_y;

          my_carry = (my_temp - my_sum) - my_y;
          my_sum   = my_temp;
      }
};

template <> 
class KahanAdder<float>
{
public:
      KahanAdder()                 : my_sum(0), my_carry(0)  {}
      KahanAdder(float e)          : my_sum(e), my_carry(0)  {}
      KahanAdder(float e, float c) : my_sum(e), my_carry(c)  {}

      // return the working sum.
      const float &sum() const { return my_sum; }

      // return the current carry.
      const float &carry() const { return my_carry; }

      operator const float&() const { return my_sum; }

      template<class S>
      KahanAdder<float> &operator =  (const S &g) 
      { my_sum = g; my_carry = 0; return *this; }

      template<class S>
      KahanAdder<float> &operator += (const S &g)  { add(g); return *this; }
      template<class S>
      KahanAdder<float> &operator -= (const S &g)  { add(-g); return *this; }
      template<class S>
      KahanAdder<float> &operator *= (const S &g)  { *this = ((float)*this) * g; return *this; }
      template<class S>
      KahanAdder<float> &operator /= (const S &g)  { *this = ((float)*this) / g; return *this; }

private:
      float my_sum;   // The current working sum.
      float my_carry; // The carry from the previous operation

      // After each add operation, the carry will contain the additional
      // bits that could be left out from the previous operation.
      void add(const float &datum) 
      {
          float my_y     = datum - my_carry;
          float my_temp  = my_sum + my_y;

          my_carry = (my_temp - my_sum) - my_y;
          my_sum   = my_temp;
      }
};

template<class T, class S>
inline KahanAdder<T> operator + (const S &c, const KahanAdder<T> &k)  
{
  KahanAdder<T> temp(k);
  temp += c;
  return temp;
}

template<class T, class S>
inline KahanAdder<T> operator + (const KahanAdder<T> &k, const S &c)
{
  KahanAdder<T> temp(k);
  temp += c;
  return temp;
}

template<class T, class S>
inline KahanAdder<T> operator + (const KahanAdder<S> &k1, const KahanAdder<T> &k2)  
{
  KahanAdder<T> temp(k1);
  temp += k2;
  return temp;
}

template<class T, class S>
inline KahanAdder<T> operator - (const S &c, const KahanAdder<T> &k)  
{
  KahanAdder<T> temp(k);
  temp -= c;
  return temp;
}

template<class T, class S>
inline KahanAdder<T> operator - (const KahanAdder<S> &k1, const KahanAdder<T> &k2)  
{
  KahanAdder<T> temp(k1);
  temp -= k2;
  return temp;
}

template<class T, class S>
inline KahanAdder<T> operator - (const KahanAdder<T> &k, const S &c)
{
  KahanAdder<T> temp(k);
  temp -= c;
  return temp;
}

#endif
