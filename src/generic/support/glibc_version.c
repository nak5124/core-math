#include <stdio.h>
#include <gnu/libc-version.h>

int main() {
  printf ("GNU libc version: %s\n", gnu_get_libc_version ());
  printf ("GNU libc release: %s\n", gnu_get_libc_release ());
  return 0;
}
