// Copyright (c) 2013 MIT License

#include "./poly.h"
#include <stdint.h>
#include "./temp_space.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
//#include <cilk/cilk_stub.h>

void karatsuba(int *c, const int *a, const int *b, size_t n);

Polynomial *polynomial_new(int size) {
  assert(size > 0);
  Polynomial *polynomial = malloc(sizeof(polynomial));
  if (polynomial == NULL) {
    return NULL;
  }
  polynomial->size = size;
  polynomial->coefficients = malloc(size * sizeof(int));
  return polynomial;
}

void polynomial_delete(Polynomial *polynomial) {
  free(polynomial->coefficients);
  free(polynomial);
}

Polynomial *multiply(Polynomial *p1, Polynomial *p2) {
  Polynomial *p3 = polynomial_new(p2->size + p1->size);
  for (int i = 0; i < p1->size; i++) {
    for (int j = 0; j < p2->size; j++) {
      p3->coefficients[i+j] += p1->coefficients[i] * p2->coefficients[j];
    }
  }
  return p3;
}

bool equals(Polynomial *p1, Polynomial *p2) {
  int i = 0;
  while (i < p1->size && i < p2->size) {
    if (p1->coefficients[i] != p2->coefficients[i]) {
      return false;
    }
    i++;
  }
  while (i < p1->size) {
    if (p1->coefficients[i] != 0) {
      return false;
    }
    i++;
  }
  while (i < p2->size) {
    if (p2->coefficients[i] != 0) {
      return false;
    }
    i++;
  }
  return true;
}

Polynomial *multiply_karat(Polynomial *p1, Polynomial *p2) {
  int n = ((p2->size > p1->size) ?
           p2->size : p1->size);
  int *a = (int *) calloc(n, sizeof(int));
  int *b = (int *) calloc(n, sizeof(int));
  int *c = (int *) calloc(2 * n, sizeof(int));

  cilk_for (size_t i = 0; i < n; i++) {
    if (i < p1->size) {
      a[i] = p1->coefficients[i];
    }
    if (i < p2->size) {
      b[i] = p2->coefficients[i];
    }
  }
  karatsuba(c, a, b, n);

  Polynomial *p3 = polynomial_new(2 * n);
  for (size_t i = 0; i < p3->size; i++) {
    p3->coefficients[i] = c[i];
  }
  free(a);
  free(b);
  free(c);

  return p3;
}

// Returns a string representation of the polynomial
void printString(Polynomial *p) {
  for (int i = p->size-1; i >= 0; i--) {
    if (p->coefficients[i] == 0) {
      continue;
    }
    printf("%d", p->coefficients[i]);
    if (i != 0) {
      printf("X^%d + ", i);
    }
  }
}

// Polynomial multiplication C=A*B that takes quadratic time.
// A = a[0:n]
// B = b[0:n]
// C = c[0:2*n-1]
// This routine is used for the "small n" basecase in Karatsuba multiplication.
void simple_mul(int *c, const int *a, const int *b, size_t n) {
  // c[0:2*n-1] = 0;
  for (size_t i = 0; i < 2 * n - 1; i++) {
    c[i] = 0;
  }
  for (size_t i = 0; i < n; i++) {
    // c[i:n] += a[i] * b[0:n];
    for (size_t j = 0; j < n; j++) {
      c[i + j] += a[i] * b[j];
    }
  }
}

// Polynomial multiplication C=A*B.
// A = a[0:n]
// B = b[0:n]
// C = c[0:2*n-1]
void karatsuba(int *c, const int *a, const int *b, size_t n) {
  if (n <= 2) {
    simple_mul(c, a, b, n);
  } else {
    size_t m = n / 2;
    // Set c[0..n-1] = $t_0$
    cilk_spawn karatsuba(c, a, b, m);

    // Set c[2*m..2*m+n-1] = $t_2$
    karatsuba(c + 2 * m, a + m, b + m, n - m);
    cilk_sync;

    // Allocate some temporary space.
    int *a_sum = temp_space(4 * (n - m));
    int *b_sum = a_sum + (n - m);
    int *t = b_sum + (n - m);

    for (int i = 0; i < m; i++) {
      a_sum[i] = a[i] + a[i + m];
      b_sum[i] = b[i] + b[i + m];
    }

    // Handle edge case, if n is odd.
    if (n & 1) {
      a_sum[m] = a[2 * m];
      b_sum[m] = b[2 * m];
    }

    // Set t = $t_1$
    karatsuba(t, a_sum, b_sum, n - m);

    // Set t = $t_1 - t_0 - t_2$
    for(int i = 0; i < 2 * m - 1; i++) {
      t[i] -= c[i] + c[2 * m + i];
    }

    // Add $(t_1 - t_0 - t_2) K$ into final product.
    c[2 * m - 1] = 0;

    // Set c[m..3*m-1] += t[0..2*m-1]
    for(int i = 0; i < 2 * m - 1; i++) {
      c[m + i] += t[i];
    }

    // Handle edge case, if n is odd.
    if (n & 1) {
      // c[3*m-1:2] += t[2*m-1:2] - c[4*m-1:2];
      c[3 * m - 1] += t[2 * m - 1] - c[4 * m - 1];
      c[3 * m] += t[2 * m] - c[4 * m];
    }

    // Clean up after ourselves.
    free(a_sum);
  }
}
