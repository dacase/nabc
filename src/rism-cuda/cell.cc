#include <stdlib.h>
#include <iostream>
#include "cell.h"
using namespace std;

#define MAX_DR 0.5

void Cell :: setup() {
  volume = box[0] * box[1] * box[2];
  ngrid = grid[0] * grid[1] * grid[2];
  dv = volume / ngrid;
  dr[0] = box[0] / grid[0];
  dr[1] = box[1] / grid[1];
  dr[2] = box[2] / grid[2];
  shift[0] = shift[1] = shift[2] = 0.0;
  if (dr[0] > MAX_DR || dr[1] > MAX_DR || dr[2] > MAX_DR) {
    cout << "##########################################" << endl;
    cout << "WARRING: Grid spacing is greater than "
	 << MAX_DR << "." << endl;
    cout << "##########################################" << endl;
  }
}
