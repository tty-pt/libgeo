/*
 * 3d_world.c - Realistic 3D voxel world example for libgeo
 *
 * Demonstrates a practical use case: simple voxel/block world
 * - Chunk-based world structure
 * - Different block types
 * - View frustum / region queries
 * - Block placement and destruction
 * - World generation
 *
 * Compile: make 3d_world
 * Run: ./3d_world
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ttypt/geo.h>
#include <ttypt/qmap.h>

// Block types
typedef enum {
	BLOCK_AIR = 0,
	BLOCK_GRASS = 1,
	BLOCK_DIRT = 2,
	BLOCK_STONE = 3,
	BLOCK_WOOD = 4,
	BLOCK_LEAVES = 5,
	BLOCK_WATER = 6,
	BLOCK_SAND = 7,
} BlockType;

const char *block_names[] = {
	"Air", "Grass", "Dirt", "Stone", "Wood", "Leaves", "Water", "Sand"
};

// World configuration
#define CHUNK_SIZE 16
#define WORLD_HEIGHT 64
#define RENDER_DISTANCE 3  // chunks

// Function prototypes
void generate_terrain(uint32_t db, int16_t chunk_x, int16_t chunk_z);
void generate_tree(uint32_t db, int16_t x, int16_t y, int16_t z);
void print_world_stats(uint32_t db);
void render_view(uint32_t db, int16_t player_x, int16_t player_y, int16_t player_z);
uint32_t count_blocks_in_region(uint32_t db, int16_t x, int16_t y, int16_t z,
                                  uint16_t w, uint16_t h, uint16_t d);

int main(void)
{
	printf("=== Libgeo 3D Voxel World Example ===\n\n");

	// Initialize
	srand(time(NULL));
	geo_init();

	// Create world database
	printf("Creating voxel world database...\n");
	uint32_t world = geo_open("world.db", "voxels", 0xFFFF);  // 65k capacity
	printf("World database created\n\n");

	// Generate a 3x3 chunk area
	printf("Generating terrain (3x3 chunks = %dx%d blocks)...\n",
	       3 * CHUNK_SIZE, 3 * CHUNK_SIZE);
	
	for (int16_t cx = -1; cx <= 1; cx++) {
		for (int16_t cz = -1; cz <= 1; cz++) {
			generate_terrain(world, cx, cz);
		}
	}
	printf("Terrain generation complete\n\n");

	// Print statistics
	print_world_stats(world);

	// Simulate player views from different positions
	printf("\n--- Simulating Player Views ---\n\n");

	printf("View 1: Player at origin (0, 10, 0)\n");
	render_view(world, 0, 10, 0);

	printf("\nView 2: Player at (16, 5, 16)\n");
	render_view(world, 16, 5, 16);

	printf("\nView 3: Player underground (-8, 2, -8)\n");
	render_view(world, -8, 2, -8);

	// Demonstrate block placement
	printf("\n--- Block Placement/Destruction ---\n\n");

	int16_t build_pos[3] = {0, 15, 0};
	printf("Building structure at (%d,%d,%d):\n",
	       build_pos[0], build_pos[1], build_pos[2]);

	// Build a small platform
	for (int16_t x = -2; x <= 2; x++) {
		for (int16_t z = -2; z <= 2; z++) {
			int16_t pos[3] = {build_pos[0] + x, build_pos[1], build_pos[2] + z};
			geo_put(world, pos, BLOCK_WOOD, 3);
		}
	}
	printf("  - Built 5x5 wooden platform\n");

	// Add walls
	for (int16_t y = 1; y <= 3; y++) {
		for (int16_t x = -2; x <= 2; x++) {
			int16_t pos1[3] = {build_pos[0] + x, build_pos[1] + y, build_pos[2] - 2};
			int16_t pos2[3] = {build_pos[0] + x, build_pos[1] + y, build_pos[2] + 2};
			geo_put(world, pos1, BLOCK_STONE, 3);
			geo_put(world, pos2, BLOCK_STONE, 3);
		}
	}
	printf("  - Built stone walls\n");

	// Count blocks in structure
	uint32_t structure_blocks = count_blocks_in_region(world,
	                                                     build_pos[0] - 2,
	                                                     build_pos[1],
	                                                     build_pos[2] - 2,
	                                                     5, 4, 5);
	printf("  - Structure contains %u blocks\n", structure_blocks);

	// Destroy some blocks
	printf("\nDestroying blocks in region (5,0,5) to (10,10,10):\n");
	int destroyed = 0;
	
	for (int16_t x = 5; x < 10; x++) {
		for (int16_t y = 0; y < 10; y++) {
			for (int16_t z = 5; z < 10; z++) {
				int16_t pos[3] = {x, y, z};
				uint32_t block = geo_get(world, pos, 3);
				if (block != GEO_MISS && block != BLOCK_AIR) {
					geo_del(world, pos, 3);
					destroyed++;
				}
			}
		}
	}
	printf("  - Destroyed %d blocks\n", destroyed);

	// Final stats
	printf("\n--- Final World Stats ---\n\n");
	print_world_stats(world);

	// Save world
	printf("\nSaving world to disk...\n");
	qmap_save();
	printf("World saved to world.db:voxels\n");

	printf("\n=== Example Complete ===\n");
	printf("\nThis example demonstrated:\n");
	printf("- Procedural terrain generation\n");
	printf("- Chunk-based world structure\n");
	printf("- View frustum queries (render distance)\n");
	printf("- Block placement and destruction\n");
	printf("- Efficient spatial queries\n");
	printf("- World persistence\n");
	printf("\nWorld file: world.db (can be deleted with: rm world.db)\n");

	qmap_close(world);
	return 0;
}

void generate_terrain(uint32_t db, int16_t chunk_x, int16_t chunk_z)
{
	int16_t base_x = chunk_x * CHUNK_SIZE;
	int16_t base_z = chunk_z * CHUNK_SIZE;

	for (int16_t lx = 0; lx < CHUNK_SIZE; lx++) {
		for (int16_t lz = 0; lz < CHUNK_SIZE; lz++) {
			int16_t x = base_x + lx;
			int16_t z = base_z + lz;

			// Simple height map (pseudo-random based on position)
			int height = 4 + ((abs(x * 7 + z * 13)) % 4);

			// Generate layers
			for (int16_t y = 0; y < height; y++) {
				int16_t pos[3] = {x, y, z};
				BlockType block;

				if (y == height - 1) {
					block = BLOCK_GRASS;  // Top layer
				} else if (y > height - 4) {
					block = BLOCK_DIRT;   // Dirt layer
				} else {
					block = BLOCK_STONE;  // Stone underground
				}

				geo_put(db, pos, block, 3);
			}

			// Occasionally place water in low areas
			if (height < 5 && (rand() % 10) < 2) {
				int16_t pos[3] = {x, height, z};
				geo_put(db, pos, BLOCK_WATER, 3);
			}

			// Occasionally place trees on grass
			if (height >= 5 && (rand() % 20) == 0) {
				generate_tree(db, x, height, z);
			}

			// Add sand near water level
			if (height == 4) {
				int16_t pos[3] = {x, height - 1, z};
				geo_put(db, pos, BLOCK_SAND, 3);
			}
		}
	}
}

void generate_tree(uint32_t db, int16_t x, int16_t y, int16_t z)
{
	// Tree trunk (3 blocks high)
	for (int16_t dy = 0; dy < 3; dy++) {
		int16_t pos[3] = {x, y + dy, z};
		geo_put(db, pos, BLOCK_WOOD, 3);
	}

	// Leaves (simple 3x3x2 cube on top)
	for (int16_t dx = -1; dx <= 1; dx++) {
		for (int16_t dz = -1; dz <= 1; dz++) {
			for (int16_t dy = 3; dy < 5; dy++) {
				if (dx == 0 && dz == 0 && dy == 3) continue;  // Skip trunk
				int16_t pos[3] = {x + dx, y + dy, z + dz};
				geo_put(db, pos, BLOCK_LEAVES, 3);
			}
		}
	}
}

void print_world_stats(uint32_t db)
{
	printf("World Statistics:\n");

	// Count blocks by type
	uint32_t type_counts[8] = {0};
	uint32_t total = 0;

	// Sample entire generated area
	int16_t start[3] = {-CHUNK_SIZE, 0, -CHUNK_SIZE};
	uint16_t lengths[3] = {CHUNK_SIZE * 3, WORLD_HEIGHT, CHUNK_SIZE * 3};
	uint32_t iter = geo_iter(db, start, lengths, 3);

	int16_t point[3];
	uint32_t value;

	while (geo_next(point, &value, iter)) {
		if (value < 8) {
			type_counts[value]++;
		}
		total++;
	}

	printf("  Total blocks: %u\n", total);
	for (int i = 1; i < 8; i++) {
		if (type_counts[i] > 0) {
			printf("  - %s: %u (%.1f%%)\n",
			       block_names[i],
			       type_counts[i],
			       100.0 * type_counts[i] / total);
		}
	}
}

void render_view(uint32_t db, int16_t player_x, int16_t player_y, int16_t player_z)
{
	// Query blocks within render distance
	int16_t view_start[3] = {
		player_x - RENDER_DISTANCE * CHUNK_SIZE,
		player_y - 10,
		player_z - RENDER_DISTANCE * CHUNK_SIZE
	};
	uint16_t view_size[3] = {
		RENDER_DISTANCE * CHUNK_SIZE * 2,
		20,
		RENDER_DISTANCE * CHUNK_SIZE * 2
	};

	uint32_t iter = geo_iter(db, view_start, view_size, 3);

	int16_t point[3];
	uint32_t value;
	uint32_t visible_blocks = 0;

	while (geo_next(point, &value, iter)) {
		visible_blocks++;
	}

	printf("  Player position: (%d, %d, %d)\n",
	       player_x, player_y, player_z);
	printf("  Render distance: %d chunks (%d blocks)\n",
	       RENDER_DISTANCE, RENDER_DISTANCE * CHUNK_SIZE);
	printf("  Visible blocks: %u\n", visible_blocks);

	// Check block at player position
	int16_t player_pos[3] = {player_x, player_y, player_z};
	uint32_t player_block = geo_get(db, player_pos, 3);

	if (player_block == GEO_MISS || player_block == BLOCK_AIR) {
		printf("  Block at player: Air (can move)\n");
	} else if (player_block < 8) {
		printf("  Block at player: %s (solid)\n", block_names[player_block]);
	}
}

uint32_t count_blocks_in_region(uint32_t db, int16_t x, int16_t y, int16_t z,
                                  uint16_t w, uint16_t h, uint16_t d)
{
	int16_t start[3] = {x, y, z};
	uint16_t lengths[3] = {w, h, d};
	uint32_t iter = geo_iter(db, start, lengths, 3);

	uint32_t count = 0;
	int16_t point[3];
	uint32_t value;

	while (geo_next(point, &value, iter)) {
		if (value != BLOCK_AIR) {
			count++;
		}
	}

	return count;
}
