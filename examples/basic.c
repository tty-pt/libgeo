/*
 * basic.c - Basic libgeo usage example
 *
 * Demonstrates fundamental operations:
 * - Initializing libgeo
 * - Opening an in-memory database
 * - Storing values at coordinates (geo_put)
 * - Retrieving values (geo_get)
 * - Deleting entries (geo_del)
 * - Handling missing entries (GEO_MISS)
 *
 * Compile: make basic
 * Run: ./basic
 */

#include <stdio.h>
#include <ttypt/geo.h>

int main(void)
{
	printf("=== Libgeo Basic Example ===\n\n");

	// Step 1: Initialize the geo subsystem
	// This MUST be called before any other geo functions
	printf("1. Initializing libgeo...\n");
	geo_init();

	// Step 2: Open an in-memory spatial database
	// NULL, NULL = in-memory only (no file persistence)
	// 0xFF = mask (256-entry initial table, auto-grows)
	printf("2. Opening in-memory database (capacity: 256)...\n");
	uint32_t db = geo_open(NULL, NULL, 0xFF);
	printf("   Database handle: %u\n\n", db);

	// Step 3: Store values at various 3D coordinates
	printf("3. Storing values at coordinates:\n");
	
	int16_t pos1[3] = {10, 20, 30};
	geo_put_3(db, pos1, 42);
	printf("   Stored value 42 at (%d, %d, %d)\n", 
	       pos1[0], pos1[1], pos1[2]);

	int16_t pos2[3] = {-5, 0, 15};
	geo_put_3(db, pos2, 100);
	printf("   Stored value 100 at (%d, %d, %d)\n",
	       pos2[0], pos2[1], pos2[2]);

	int16_t pos3[3] = {0, 0, 0};
	geo_put_3(db, pos3, 1);
	printf("   Stored value 1 at origin (0, 0, 0)\n\n");

	// Step 4: Retrieve values from coordinates
	printf("4. Retrieving values:\n");
	
	uint32_t val = geo_get_3(db, pos1);
	printf("   Value at (%d, %d, %d) = %u\n",
	       pos1[0], pos1[1], pos1[2], val);

	val = geo_get_3(db, pos2);
	printf("   Value at (%d, %d, %d) = %u\n",
	       pos2[0], pos2[1], pos2[2], val);

	val = geo_get_3(db, pos3);
	printf("   Value at (%d, %d, %d) = %u\n\n",
	       pos3[0], pos3[1], pos3[2], val);

	// Step 5: Query a coordinate with no stored value
	printf("5. Querying empty coordinate:\n");
	int16_t empty[3] = {100, 100, 100};
	val = geo_get_3(db, empty);
	
	if (val == GEO_MISS) {
		printf("   No value at (%d, %d, %d) - got GEO_MISS\n\n",
		       empty[0], empty[1], empty[2]);
	} else {
		printf("   Unexpected: found value %u\n\n", val);
	}

	// Step 6: Update an existing value
	printf("6. Updating existing value:\n");
	printf("   Old value at (%d, %d, %d) = %u\n",
	       pos1[0], pos1[1], pos1[2], geo_get_3(db, pos1));
	
	geo_set_3(db, pos1, 99);  // Replace 42 with 99
	printf("   New value at (%d, %d, %d) = %u\n\n",
	       pos1[0], pos1[1], pos1[2], geo_get_3(db, pos1));

	// Step 7: Delete an entry
	printf("7. Deleting entry:\n");
	geo_del_3(db, pos2);
	printf("   Deleted entry at (%d, %d, %d)\n",
	       pos2[0], pos2[1], pos2[2]);
	
	val = geo_get_3(db, pos2);
	if (val == GEO_MISS) {
		printf("   Confirmed: coordinate is now empty\n\n");
	}

	// Step 8: Demonstrate 2D usage
	printf("8. Using 2D coordinates:\n");
	int16_t pos2d[2] = {50, 60};
	geo_put_2(db, pos2d, 777);  // 2D point
	printf("   Stored value 777 at (%d, %d) in 2D\n",
	       pos2d[0], pos2d[1]);
	
	val = geo_get_2(db, pos2d);  // 2D point
	printf("   Retrieved value %u from (%d, %d)\n\n",
	       val, pos2d[0], pos2d[1]);

	// Note: In a real application with file-backed databases,
	// you would call qmap_save() to persist data mid-execution,
	// or rely on automatic save at process exit.

	printf("=== Example Complete ===\n");
	printf("Database remains in memory until process exit.\n");
	printf("File-backed databases auto-save on exit.\n");

	return 0;
}
