// Parallel Fast Multipole Method
// Copyright (C) 2001-2026  Victor Anisimov
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#define _GNU_SOURCE 1
#include <stdio.h>
#include <sched.h>
#include <stdlib.h>
#include <numa.h>
#include <numaif.h>
#include <errno.h>
//#include <sycl/sycl.hpp>

#define GB (1024*1024*1024)
#define NODEMAX 4

// Return free GPU memory
//long long stat_gpumem() {
//  sycl::queue Q;
//  return Q.get_device().get_info<sycl::ext::intel::info::device::free_memory>();
//}


int getcpuid() {
  unsigned int cpuid = 0, socketid = 0;
  int int_cpuid;
  getcpu( &cpuid, &socketid );
  int_cpuid = (int) cpuid;
  return int_cpuid;
}


void* allocate_hbm_memory(int64_t size) {
  unsigned int cpuid = 0, socketid = 0;
  long long node_size[NODEMAX], free_size[NODEMAX];
  int my_hbm_numa_node;

  // ALCF Aurora
  getcpu( &cpuid, &socketid );
  if(socketid == 0)
    my_hbm_numa_node = 2;  // CPU socket 0
  else
    my_hbm_numa_node = 3;  // CPU socket 1

  // Get the amount of free memory
  for (int i = 0; i <= NODEMAX - 1; i++) {
    node_size[i] = numa_node_size( i, &(free_size[i]) );
    if(node_size[i] == -1) fprintf(stderr, "Error getting size for node %d\n", i);
  }

  // check if there is sufficient free memory
  if((long long) free_size[my_hbm_numa_node] - size > (long long) 0) {
    printf("Allocated %5.2f GB of HBM memory\n", (float)size/GB);
    printf("CPU %3d, Socket %d, HBM nodeid %d, Free HBM memory: %5.2f GB\n", 
      cpuid, socketid, my_hbm_numa_node, (float) (free_size[my_hbm_numa_node]-size)/GB);
    return numa_alloc_onnode(size, my_hbm_numa_node);
  }
  else {
    printf("Error: unable to allocate %5.2f GB of HBM memory\n", (float)size/GB);
    printf("CPU %3d, Socket %d, HBM nodeid %d, Free HBM memory: %5.2f GB\n", 
      cpuid, socketid, my_hbm_numa_node, (float) free_size[my_hbm_numa_node]/GB);
    return NULL;
  }
}


void free_hbm_memory(void** hbm_buffer, int64_t size) {
  numa_free(*hbm_buffer, size);
  *hbm_buffer = NULL;
  return;
}


void stat_numa_memory() {
  unsigned int cpuid = 0, socketid = 0;
  long long node_size[NODEMAX], free_size[NODEMAX];
  int my_hbm_numa_node;

  // ALCF Aurora
  getcpu( &cpuid, &socketid );
  if(socketid == 0)
    my_hbm_numa_node = 2;  // CPU socket 0
  else
    my_hbm_numa_node = 3;  // CPU socket 1

  // Get the amount of free memory
  for (int i = 0; i <= NODEMAX - 1; i++) {
    node_size[i] = numa_node_size( i, &(free_size[i]) );
    if(node_size[i] == -1) fprintf(stderr, "Error getting size for node %d\n", i);
  }

  printf("CPU %3d, SOCKET %d:  Free memory: %12.2lf    %12.2lf    %12.2lf    %12.2lf\n",
    cpuid, socketid,
    (float)free_size[0]/GB, (float)free_size[1]/GB,
    (float)free_size[2]/GB, (float)free_size[3]/GB );

  return;
}


long long stat_memory(int numa_node) {
  long long node_mem_size, free_mem_size;  // in bytes
  node_mem_size = numa_node_size( numa_node, &free_mem_size);
  return free_mem_size;
}


long long stat_allmemory() {
  long long node_mem_size, free_allmem_size, free_nodemem_size;  // in bytes

  // Get the amount of free memory
  free_allmem_size = 0;
  for (int numa_node = 0; numa_node <= NODEMAX - 1; numa_node++) {
    node_mem_size = numa_node_size( numa_node, &free_nodemem_size);
    free_allmem_size += free_nodemem_size;
  }
  return free_allmem_size;
}

