/*
 * iteration.c - Spatial iteration example for libgeo
 *
 * Demonstrates spatial queries and iteration:
 * - Populating a spatial database with multiple points
 * - Using geo_iter() to query rectangular regions
 * - Using geo_next() to iterate through results
 * - Observing Morton code ordering effects
 * - Querying empty regions
 *
 * Compile: make iteration
 * Run: ./iteration
 */

#include <stdio.h>
#include <ttypt/geo.h>

int main(void)
{
	printf("=== Libgeo Iteration Example ===\n\n");

	// Initialize and create database
	geo_init();
	uint32_t db = geo_open(NULL, NULL, 0xFFF);  // 4096 capacity
	printf("Created spatial database\n\n");

	// Step 1: Populate a 10x10x10 grid with values
	printf("1. Populating 10x10x10 grid (storing coordinate sum):\n");
	int count = 0;
	
	for (int16_t x = 0; x < 10; x++) {
		for (int16_t y = 0; y < 10; y++) {
			for (int16_t z = 0; z < 10; z++) {
				int16_t pos[3] = {x, y, z};
				uint32_t value = x + y + z;  // Sum of coordinates
				geo_put(db, pos, value, 3);
				count++;
			}
		}
	}
	printf("   Stored %d points\n\n", count);

	// Step 2: Query a small region
	printf("2. Querying region [0,0,0] to [3,3,3] (3x3x3 cube):\n");
	
	int16_t start1[3] = {0, 0, 0};
	uint16_t lengths1[3] = {3, 3, 3};  // 3x3x3 = 27 points
	uint32_t iter1 = geo_iter(db, start1, lengths1, 3);
	
	int16_t point[3];
	uint32_t value;
	int result_count = 0;
	
	printf("   Results (Morton/Z-order):\n");
	while (geo_next(point, &value, iter1)) {
		printf("   (%2d,%2d,%2d) = %2u", 
		       point[0], point[1], point[2], value);
		
		// Show first few to demonstrate ordering
		if (result_count < 10) {
			printf("\n");
		} else if (result_count == 10) {
			printf("   ... (showing first 10)\n");
		}
		result_count++;
	}
	printf("   Total: %d points found\n\n", result_count);

	// Step 3: Query a slice (2D plane)
	printf("3. Querying Z=5 plane, region [2,2] to [8,8]:\n");
	
	int16_t start2[3] = {2, 2, 5};
	uint16_t lengths2[3] = {6, 6, 1};  // 6x6x1 = 36 points
	uint32_t iter2 = geo_iter(db, start2, lengths2, 3);
	
	result_count = 0;
	uint32_t sum = 0;
	
	while (geo_next(point, &value, iter2)) {
		sum += value;
		result_count++;
	}
	printf("   Found %d points, sum of values = %u\n", result_count, sum);
	printf("   Average value = %.1f\n\n", (float)sum / result_count);

	// Step 4: Query overlapping region
	printf("4. Querying region [5,5,5] to [15,15,15]:\n");
	printf("   (extends beyond our 10x10x10 grid)\n");
	
	int16_t start3[3] = {5, 5, 5};
	uint16_t lengths3[3] = {10, 10, 10};  // 10x10x10 potential
	uint32_t iter3 = geo_iter(db, start3, lengths3, 3);
	
	result_count = 0;
	while (geo_next(point, &value, iter3)) {
		result_count++;
	}
	printf("   Found %d points (only within grid bounds)\n", result_count);
	printf("   Expected: 5x5x5 = 125 points\n\n");

	// Step 5: Query empty region
	printf("5. Querying empty region [100,100,100] to [110,110,110]:\n");
	
	int16_t start4[3] = {100, 100, 100};
	uint16_t lengths4[3] = {10, 10, 10};
	uint32_t iter4 = geo_iter(db, start4, lengths4, 3);
	
	result_count = 0;
	while (geo_next(point, &value, iter4)) {
		result_count++;
	}
	printf("   Found %d points (empty region)\n\n", result_count);

	// Step 6: Single point query using iteration
	printf("6. Using iteration to check specific point [7,8,9]:\n");
	
	int16_t start5[3] = {7, 8, 9};
	uint16_t lengths5[3] = {1, 1, 1};  // Single point
	uint32_t iter5 = geo_iter(db, start5, lengths5, 3);
	
	if (geo_next(point, &value, iter5)) {
		printf("   Point (%d,%d,%d) exists with value %u\n",
		       point[0], point[1], point[2], value);
		printf("   (Expected: 7+8+9 = 24)\n\n");
	} else {
		printf("   Point not found\n\n");
	}

	// Step 7: Demonstrate 2D iteration
	printf("7. 2D example - populate and query:\n");
	
	// Clear and repopulate with 2D data
	uint32_t db2d = geo_open(NULL, NULL, 0xFF);
	
	// Create a checkerboard pattern
	for (int16_t x = 0; x < 8; x++) {
		for (int16_t y = 0; y < 8; y++) {
			if ((x + y) % 2 == 0) {  // Checkerboard
				int16_t pos[2] = {x, y};
				geo_put(db2d, pos, 1, 2);  // dim=2
			}
		}
	}
	
	// Query the 2D region
	int16_t start2d[2] = {0, 0};
	uint16_t lengths2d[2] = {8, 8};
	uint32_t iter2d = geo_iter(db2d, start2d, lengths2d, 2);
	
	result_count = 0;
	while (geo_next(point, &value, iter2d)) {
		result_count++;
	}
	printf("   Checkerboard: stored %d points in 8x8 grid\n", result_count);
	printf("   (Expected: 32 squares on checkerboard)\n\n");

	// Step 8: Demonstrate sparse data
	printf("8. Sparse data example:\n");
	uint32_t db_sparse = geo_open(NULL, NULL, 0xFF);
	
	// Store only a few points in large region
	int16_t sparse_points[][3] = {
		{-1000, -1000, -1000},
		{0, 0, 0},
		{1000, 1000, 1000},
	};
	
	for (int i = 0; i < 3; i++) {
		geo_put(db_sparse, sparse_points[i], i + 1, 3);
	}
	
	// Query large region
	int16_t start_sparse[3] = {-2000, -2000, -2000};
	uint16_t lengths_sparse[3] = {4000, 4000, 4000};
	uint32_t iter_sparse = geo_iter(db_sparse, start_sparse, lengths_sparse, 3);
	
	printf("   Query region: 4000x4000x4000 (64 billion potential points)\n");
	result_count = 0;
	while (geo_next(point, &value, iter_sparse)) {
		printf("   Found: (%d,%d,%d) = %u\n",
		       point[0], point[1], point[2], value);
		result_count++;
	}
	printf("   Total found: %d (sparse database efficient!)\n\n", result_count);

	printf("=== Example Complete ===\n");
	printf("\nKey observations:\n");
	printf("- Iteration returns results in Morton/Z-order (not spatial order)\n");
	printf("- Empty cells are automatically skipped\n");
	printf("- Regions beyond stored data return no results\n");
	printf("- Sparse data is handled efficiently\n");
	printf("- Works with any dimensionality (2D, 3D, etc.)\n");

	return 0;
}
