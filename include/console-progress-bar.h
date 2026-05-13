/*
  Console Progress Bar: A C code for printing the progress of numeric
  iterations on the console

  Repository: <https://github.com/jodesarro/console-progress-bar>
  License: Refer to the LICENSE file in the Repository
  Language standard: C99
  Last updated: 2026-05-13

  Description: Code for handling progress.
*/
#ifndef CONSOLE_PROGRESS_BAR_H
#define CONSOLE_PROGRESS_BAR_H

#include <stdio.h>

/*
  Print a progress bar for an iteration process.

  Parameters:
  - iteration_current, the number of the current iteration.
  - iteration_total, the total number of iterations.

  Implementation: It prints a progress bar with a percentage evaluated through
  100 * iteration_current / iteration_total.
*/
static inline void print_progress_bar(long int iteration_current,
                                      long int iteration_total) {
  if (iteration_current == 0) {
    printf("\rProgress: [          ]   0%%");
    fflush(stdout);
  } else if (iteration_current == iteration_total) {
    printf("\rProgress: [==========] 100%%");
    fflush(stdout);
  } else {
    int percent = (int)((100.01 * iteration_current) / iteration_total);
    int pos = (int)((10.01 * iteration_current) / iteration_total);
    printf("\rProgress: [");
    for (int i = 0; i < 10; i++) {
      if (i < pos)
        printf("=");
      else
        printf(" ");
    }
    printf("] %3d%%", percent);
    fflush(stdout);
  }
}

/*
  Print a progress bar at a fixed percentage interval.

  Parameters:
  - iteration_current, the number of the current iteration.
  - iteration_total, the total number of iterations.
  - percent, the percentage interval.

  Implementation: It calls print_progress_bar(iteration_current,
  iteration_total) at every percent of percentage.
*/
static inline void print_progress_bar_every_percent(long int iteration_current,
                                                    long int iteration_total,
                                                    int percent) {
  if (iteration_current % (iteration_total * ((long int)percent) / 100L) ==
      0L) {
    print_progress_bar(iteration_current, iteration_total);
  }
}

/*
  Print an empty progress bar for a 0% progress.

  Implementation: It calls print_progress_bar(0, 0).
*/
static inline void print_progress_bar_empty() { print_progress_bar(0, 0); }

/*
  Print a full progress bar for a 100% progress followed by a line break.

  Implementation: It calls print_progress_bar(1, 1) and then prints a line
  break.
*/
static inline void print_progress_bar_full() {
  print_progress_bar(1, 1);
  printf("\n");
}

#endif /* CONSOLE_PROGRESS_BAR_H */