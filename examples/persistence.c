/*
 * persistence.c - File persistence example for libgeo
 *
 * Demonstrates file-backed databases:
 * - Opening databases with file persistence
 * - Storing and retrieving data
 * - Automatic persistence at exit
 * - Reopening to verify data was saved
 * - Multiple databases in a single file
 *
 * Compile: make persistence
 * Run: ./persistence
 */

#include <stdio.h>
#include <unistd.h>
#include <ttypt/geo.h>
#include <ttypt/qmap.h>

void example_initial_save(void);
void example_verify_load(void);
void example_multiple_databases(void);

int main(void)
{
	printf("=== Libgeo Persistence Example ===\n\n");

	printf("Part 1: Creating and saving data\n");
	printf("================================\n");
	example_initial_save();

	printf("\n\nPart 2: Verifying data persisted\n");
	printf("=================================\n");
	example_verify_load();

	printf("\n\nPart 3: Multiple databases in one file\n");
	printf("=======================================\n");
	example_multiple_databases();

	printf("\n=== Example Complete ===\n");
	printf("\nGenerated files:\n");
	printf("  - example_world.db\n");
	printf("  - example_multi.db\n");
	printf("\nThese files contain the persistent spatial databases.\n");
	printf("You can delete them with: rm example_*.db\n");

	return 0;
}

void example_initial_save(void)
{
	// Initialize
	geo_init();

	// Open a file-backed database
	printf("Opening file-backed database: example_world.db:main\n");
	uint32_t db = geo_open("example_world.db", "main", 0xFF);

	// Store some data representing a simple voxel world
	printf("Storing voxel data:\n");

	// Ground layer (y=0, z=0..9, x=0..9)
	for (int16_t x = 0; x < 10; x++) {
		for (int16_t z = 0; z < 10; z++) {
			int16_t pos[3] = {x, 0, z};
			geo_put(db, pos, 1, 3);  // Block type 1 = grass
		}
	}
	printf("  - Added 10x10 ground layer (type 1 = grass)\n");

	// Add some trees at specific locations
	int16_t tree_positions[][3] = {
		{3, 1, 3}, {3, 2, 3}, {3, 3, 3},  // Tree 1
		{7, 1, 7}, {7, 2, 7}, {7, 3, 7},  // Tree 2
	};
	
	for (unsigned i = 0; i < sizeof(tree_positions) / sizeof(tree_positions[0]); i++) {
		geo_put(db, tree_positions[i], 2, 3);  // Type 2 = wood
	}
	printf("  - Added 2 trees (type 2 = wood)\n");

	// Add player spawn point
	int16_t spawn[3] = {5, 1, 5};
	geo_put(db, spawn, 100, 3);  // Type 100 = spawn
	printf("  - Added player spawn at (%d,%d,%d)\n",
	       spawn[0], spawn[1], spawn[2]);

	// Explicitly save (optional - auto-saves on exit anyway)
	printf("\nCalling qmap_save() to persist data...\n");
	qmap_save();
	printf("Data saved to example_world.db\n");

	// Close database (optional - auto-closes on exit)
	qmap_close(db);
	printf("Database closed\n");
}

void example_verify_load(void)
{
	// Initialize again (simulating a new process)
	geo_init();

	// Reopen the same file - data should auto-load
	printf("Reopening example_world.db:main\n");
	uint32_t db = geo_open("example_world.db", "main", 0xFF);
	printf("Database opened, data auto-loaded from file\n\n");

	// Verify ground layer
	printf("Verifying data:\n");
	int ground_count = 0;
	for (int16_t x = 0; x < 10; x++) {
		for (int16_t z = 0; z < 10; z++) {
			int16_t pos[3] = {x, 0, z};
			uint32_t value = geo_get(db, pos, 3);
			if (value == 1) ground_count++;
		}
	}
	printf("  - Ground layer: %d/100 blocks found\n", ground_count);

	// Verify spawn point
	int16_t spawn[3] = {5, 1, 5};
	uint32_t spawn_value = geo_get(db, spawn, 3);
	if (spawn_value == 100) {
		printf("  - Spawn point: FOUND at (%d,%d,%d)\n",
		       spawn[0], spawn[1], spawn[2]);
	} else {
		printf("  - Spawn point: NOT FOUND (got %u)\n", spawn_value);
	}

	// Count all blocks using iteration
	int16_t start[3] = {-10, -10, -10};
	uint16_t lengths[3] = {30, 30, 30};
	uint32_t iter = geo_iter(db, start, lengths, 3);
	
	int total_blocks = 0;
	int16_t point[3];
	uint32_t value;
	
	while (geo_next(point, &value, iter)) {
		total_blocks++;
	}
	
	printf("  - Total blocks: %d\n", total_blocks);
	printf("\n✓ Data successfully persisted and loaded!\n");

	qmap_close(db);
}

void example_multiple_databases(void)
{
	geo_init();

	// Create multiple logical databases in the same file
	printf("Creating three databases in example_multi.db:\n");

	// Database 1: Player data
	uint32_t db_players = geo_open("example_multi.db", "players", 0xFF);
	int16_t player1[3] = {100, 50, 200};
	int16_t player2[3] = {-50, 60, -30};
	geo_put(db_players, player1, 1001, 3);  // Player ID 1001
	geo_put(db_players, player2, 1002, 3);  // Player ID 1002
	printf("  1. 'players' database: 2 player positions\n");

	// Database 2: Chunks
	uint32_t db_chunks = geo_open("example_multi.db", "chunks", 0xFF);
	for (int16_t cx = 0; cx < 5; cx++) {
		for (int16_t cz = 0; cz < 5; cz++) {
			int16_t chunk_pos[2] = {cx, cz};
			uint32_t chunk_id = cx * 1000 + cz;
			geo_put(db_chunks, chunk_pos, chunk_id, 2);  // 2D chunks
		}
	}
	printf("  2. 'chunks' database: 25 chunk locations\n");

	// Database 3: Items
	uint32_t db_items = geo_open("example_multi.db", "items", 0xFF);
	int16_t item_positions[][3] = {
		{10, 5, 20},
		{15, 3, 18},
		{-5, 10, 25},
	};
	for (unsigned i = 0; i < sizeof(item_positions) / sizeof(item_positions[0]); i++) {
		geo_put(db_items, item_positions[i], 500 + i, 3);  // Item IDs
	}
	printf("  3. 'items' database: 3 dropped items\n");

	// Save all databases
	qmap_save();
	printf("\nAll databases saved to example_multi.db\n");

	// Close
	qmap_close(db_players);
	qmap_close(db_chunks);
	qmap_close(db_items);

	// Reopen and verify
	printf("\nReopening and verifying:\n");
	
	db_players = geo_open("example_multi.db", "players", 0xFF);
	db_chunks = geo_open("example_multi.db", "chunks", 0xFF);
	db_items = geo_open("example_multi.db", "items", 0xFF);

	uint32_t p1 = geo_get(db_players, player1, 3);
	uint32_t p2 = geo_get(db_players, player2, 3);
	printf("  - Player 1: ID %u at (%d,%d,%d)\n", p1,
	       player1[0], player1[1], player1[2]);
	printf("  - Player 2: ID %u at (%d,%d,%d)\n", p2,
	       player2[0], player2[1], player2[2]);

	int16_t chunk_test[2] = {2, 3};
	uint32_t chunk = geo_get(db_chunks, chunk_test, 2);
	printf("  - Chunk at (%d,%d): ID %u\n",
	       chunk_test[0], chunk_test[1], chunk);

	uint32_t item = geo_get(db_items, item_positions[0], 3);
	printf("  - First item: ID %u at (%d,%d,%d)\n", item,
	       item_positions[0][0], item_positions[0][1], item_positions[0][2]);

	printf("\n✓ Multiple databases successfully stored in one file!\n");

	qmap_close(db_players);
	qmap_close(db_chunks);
	qmap_close(db_items);
}
