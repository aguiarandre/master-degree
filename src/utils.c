#include "native.h"

#include <stddef.h>

/* get the prime number after offset. */
uint utils_next_prime(uint offset)
{
  int prime[21];
  uint x, i, q, r;

  prime[0] = 0;
  prime[1] = 3;
  prime[2] = 5;
  prime[3] = 7;
  prime[4] = 11;
  prime[5] = 13;
  prime[6] = 17;
  prime[7] = 19;
  prime[8] = 23;
  prime[9] = 29;
  prime[10] = 31;
  prime[11] = 37;
  prime[12] = 41;
  prime[13] = 43;
  prime[14] = 47;
  prime[15] = 53;
  prime[16] = 59;
  prime[17] = 61;
  prime[18] = 67;
  prime[19] = 71;
  prime[20] = 73;

  x = 1 + 2*(7*offset/10);

  //if ( x > prime[20] ) goto label_30;
  if ( x >= prime[20] || x < prime[20] ) goto label_30;
  i = 21;
  label_10:
    if ( x >= prime[i] ) goto label_20;
    i -= 1;
    goto label_10;
  label_20:
    if ( x > prime[i] ) i += 1;
    return prime[i];
  label_30:
    for ( i = 1; i <= 20; i++ ) {
      q = x/prime[i];
      r = x - q*prime[i];
      if ( r != 0 ) goto label_40;
      x += 2;
      goto label_30;
    label_40:
      if ( q <= prime[i] ) goto label_60;
    }
  label_60:
    return x;
}


/* function for unsigned int comparison, to be used
   with qsort(). */
int compare_uint(const void* val1, const void* val2 )
{
  uint a = *(uint *)val1;
  uint b = *(uint *)val2;
  if (a > b)
      return 1;
  else if (a < b)
      return -1;
  else
      return 0;
}

/* function to remove duplicate entries in 
   a sorted uint array. Similar to c++ std::unique.
   The array is changed in place and the function
   returns the array new size. */
uint unique (uint* first, uint* last)
{
  if (first==last) return *last;

  uint size = 0;
  uint* result = first;
  while (++first != last)
  {
    if (!(*result == *first)) {
      *(++result)=*first;
      ++size;
    }
  }
  //return ++result;
  return ++size;
}

/* function to binary search an ordered array.
 * Used to look for renumbered vertices
 * indeces for GUI processing. */
uint binary_search( uint* array, uint size, uint value )
{
  uint index = -1;

  uint first = 0;
  uint last = size - 1;
  uint middle = (first+last)/2;

  while( first <= last )
  {
    if ( array[middle] < value )
        first = middle + 1;
    else if ( array[middle] == value )
    {
      index = middle;
      break;
    }
    else {
      last = middle - 1;
    }
    middle = (first + last)/2;
  }

  return index;
}
