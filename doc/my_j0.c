/*
  $ gcc -fPIC -c my_j0.c -I$MPFI/include
  $ gcc -shared -o my_j0 my_j0.o -lgmp -lmpfr
  In Sollya:
  > my_j0=library("./my_j0");
  > evaluate(my_j0(x), 2);
  0.22389077914123566805182745464994862582515448221861
  > evaluate(diff(my_j0(x)),2);
  -0.57672480775687338720244824226913708692030268971967
  > evaluate(diff(diff(my_j0(x))),2);
  6.447162473720102554939666648461991763499686264123e-2
  > evaluate(diff(diff(diff(my_j0(x)))),2);
  [-infty;infty]
 */

#include <mpfr.h>
#include "mpfi.h"

#define left(x) (&(x->left))
#define right(x) (&(x->right))

/* assumes j0 and its first and second derivatives are monotonic on op=[a,b] */
int
my_j0 (mpfi_t rop, mpfi_t op, int n)
{
  if (n == 0)
    {
      /* since we assumed j0 monotonic, it is included in (j0(a),j0(b)) */
      mpfr_j0 (left (rop), left (op), MPFR_RNDU);
      mpfr_j0 (right (rop), right (op), MPFR_RNDU);
      if (mpfr_cmp (left (rop), right (rop)) <= 0)
        {
          /* j0 is increasing, we should round j0(a) down */
          mpfr_nextbelow (left (rop));
        }
      else
        {
          /* j0 is decreasing, we should round j0(b) down and swap */
          mpfr_nextbelow (right (rop));
          mpfr_swap (left (rop), right (rop));
        }
      return 1;
    }

   if (n == 1)
     {
       /* diff(j0) = -j1, we assume j1 is monotonic on op = [a,b] */
       mpfr_j1 (left (rop), left(op), MPFR_RNDD);
       mpfr_neg (left (rop), left (rop), MPFR_RNDU);
       mpfr_j1 (right (rop), right(op), MPFR_RNDD);
       mpfr_neg (right (rop), right (rop), MPFR_RNDU);
      if (mpfr_cmp (left (rop), right (rop)) <= 0)
        {
          /* -j1 is increasing, we should round -j1(a) down */
          mpfr_nextbelow (left (rop));
        }
      else
        {
          /* -j1 is decreasing, we should round -j1(b) down and swap */
          mpfr_nextbelow (right (rop));
          mpfr_swap (left (rop), right (rop));
        }
       return 1;
   }

   if (n == 2)
     {
       /* diff(diff(j0)) = -diff(j1) = -1/2 (j0 - j2) = 1/2 (j2 - j0)
          we assume both j0, j2 and diff(diff(j0)) are monotonic on op = [a,b] */
       mpfr_t s, t, u, v;
       mpfr_init2 (s, mpfi_get_prec (rop));
       mpfr_init2 (t, mpfi_get_prec (rop));
       mpfr_init2 (u, mpfi_get_prec (rop));
       mpfr_init2 (v, mpfi_get_prec (rop));
       mpfr_j0 (s, left (op), MPFR_RNDD);      /* s <= j0(a) */
       mpfr_jn (t, 2, left (op), MPFR_RNDU);   /* t >= j2(a) */
       mpfr_j0 (u, right (op), MPFR_RNDD);     /* u <= j0(b) */
       mpfr_jn (v, 2, right (op), MPFR_RNDU);  /* v >= j2(b) */
       mpfr_sub (left (rop), t, s, MPFR_RNDU); /* left(rop) >= j2(a)-j0(a) */
       mpfr_sub (right (rop), v, u, MPFR_RNDU); /* right(rop) >= j2(b)-j0(b) */
       if (mpfr_cmp (left (rop), right (rop)) <= 0)
         {
           /* diff(diff(j0)) is increasing on [a,b], thus we round down in a */
           mpfr_nextabove (s); /* now s >= j0(a) */
           mpfr_nextbelow (t); /* now t <= j2(a) */
           mpfr_sub (left (rop), t, s, MPFR_RNDD); /* now left(rop) <= j2(a)-j0(a) */
         }
       else
         {
           /* diff(diff(j0)) is decreasing on [a,b], thus we round down in b */
           mpfr_nextabove (u); /* now u >= j0(b) */
           mpfr_nextbelow (v); /* now v <= j2(b) */
           mpfr_sub (right (rop), v, u, MPFR_RNDD); /* right(rop) <= j2(b)-j0(b) */
           mpfr_swap (left (rop), right (rop));
         }
       mpfr_div_2ui (left (rop), left (rop), 1, MPFR_RNDD); /* exact */
       mpfr_div_2ui (right (rop), right (rop), 1, MPFR_RNDU); /* exact */
       mpfr_clear (s);
       mpfr_clear (t);
       mpfr_clear (u);
       mpfr_clear (v);
       return 1;
     }

   /* larger values of n should not be useful */
   mpfr_set_inf (left (rop), -1);
   mpfr_set_inf (right (rop), 1);
   return 1;
}
