#include "cache.h"

uint32_t partitioning[4][2] = { {7,9},{8,8},{12,4}, {1,15} };  //stores the partitioning of the set for core0 and core1
uint32_t count[2]; //stores the number of blocks of each core for a particular set
uint32_t pidx = 1; // stores the index of the partition
uint32_t vic_cpu = 0; //stores the core number whose block will be evicted

static int counter[16] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
static int as = 0;
int allocations[NUM_CPUS];
int blocks_req[NUM_CPUS];
int alloc;
//int balance=NUM_WAY;
int max_mu[NUM_CPUS];
int get_max_mu(int i, int alloc, int b);
int get_max_k(int i, int alloc, int b);
int get_mu_value(int p, int a, int b);
int balance = LLC_WAY;
static int n = 0;
int arr[32][LLC_WAY];






uint32_t CACHE::find_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK* current_set, uint64_t ip, uint64_t full_addr, uint32_t type)
{
    // baseline LRU replacement policy for other caches 
    return lru_victim(cpu, instr_id, set, current_set, ip, full_addr, type);
}

void CACHE::update_replacement_state(uint32_t cpu, uint32_t set, uint32_t way, uint64_t full_addr, uint64_t ip, uint64_t victim_addr, uint32_t type, uint8_t hit)
{
    if (type == WRITEBACK) {
        if (hit) // writeback hit does not update LRU state
            return;
    }

    return lru_update(set, way);
}

//function to check number of misses
void CACHE::partition(int i, int j)
{

    if (as != 32) {
        if (i == 0 || (i - 1) % 32 == 0) {
            counter[j]++;
        }
    }
}
//Lookahead algorithm for ATD 
int* algo()
{


    for (int i = 0;i < NUM_CPUS;i++) {
        allocations[i] = 0;
    }

    //int b=NUM_WAY;
    while (balance > 0)
    {
        for (int k = 0;k < NUM_CPUS;k++) {
            alloc = allocations[k];
            max_mu[k] = get_max_mu(k, alloc, balance);
            blocks_req[k] = get_max_k(k, alloc, balance);
        }


        int ma = max_mu[0];
        int win;
        for (int m = 0;m < NUM_CPUS;m++)
        {
            if (max_mu[m] > ma) {
                ma = max_mu[m];win = m;
            }
        }

        allocations[win] += blocks_req[win];
        balance -= blocks_req[win];
    }

    return allocations;

}
int get_max_mu(int i, int alloc, int b) {  //Function to get maximum utility
    int max = 0;
    int mu;
    for (int k = 0;k <= b;k++)
    {
        mu = get_mu_value(i, alloc, alloc + k);
        if (mu > max) {
            max = mu;
        }
    }
    return max;
}
int get_max_k(int i, int alloc, int b) {  //Function to get blocks corresponding maximum utility
    int max = 0;
    int mu;
    int r;
    for (int k = 0;k <= b;k++)
    {
        mu = get_mu_value(i, alloc, alloc + k);
        if (mu > max) {
            max = mu;r = k;
        }
    }
    return r;
}
int get_mu_value(int p, int a, int b)  //Function to get utility
{
    int u = counter[b - 1] - counter[a - 1];
    return u / (b - a);
}





///////////////////////////////////////////////////////////////////////////////////////


uint32_t CACHE::lru_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK* current_set, uint64_t ip, uint64_t full_addr, uint32_t type)
{
    uint32_t way = 0;


    count[0] = 0;
    count[1] = 0;
    // counts the number of blocks of each core in the set
    for (way = 0; way < NUM_WAY; way++) {
        if (block[set][way].valid == true) {
            count[block[set][way].cpu]++;
        }
    }

    // fill invalid line first
    for (way = 0; way < NUM_WAY; way++) {
        if (block[set][way].valid == false) {

            DP(if (warmup_complete[cpu]) {
                cout << "[" << NAME << "] " << func << " instr_id: " << instr_id << " invalid set: " << set << " way: " << way;
                cout << hex << " address: " << (full_addr >> LOG2_BLOCK_SIZE) << " victim address: " << block[set][way].address << " data: " << block[set][way].data;
                cout << dec << " lru: " << block[set][way].lru << endl;
            });

            break;
        }
    }

    // LRU victim
    if (way == NUM_WAY) {
        for (way = 0; way < NUM_WAY; way++) {
            if (block[set][way].lru == NUM_WAY - 1) {

                DP(if (warmup_complete[cpu]) {
                    cout << "[" << NAME << "] " << func << " instr_id: " << instr_id << " replace set: " << set << " way: " << way;
                    cout << hex << " address: " << (full_addr >> LOG2_BLOCK_SIZE) << " victim address: " << block[set][way].address << " data: " << block[set][way].data;
                    cout << dec << " lru: " << block[set][way].lru << endl;
                });

                break;
            }
        }
    }
    // }
    if (way == NUM_WAY) {
        cerr << "[" << NAME << "] " << " no victim! set: " << set << endl;///nnnnnnnnn
        assert(0);
    }
    return way;
}

uint32_t CACHE::llc_lru_victim(uint32_t cpu, uint64_t instr_id, uint32_t set, const BLOCK* current_set, uint64_t ip, uint64_t full_addr, uint32_t type)
{
    n++;
    if (n == 32)
    {
        int* at;
        at = algo();
        partitioning[pidx][0] = at[0];
        partitioning[pidx][1] = LLC_WAY - at[0];
    }
    uint32_t way = 0;
    count[0] = 0;   // number of block of core0
    count[1] = 0;   // number of block of core1

    // counts the number of blocks of each core in the set
    for (way = 0; way < NUM_WAY; way++) {
        if (block[set][way].valid == true) {
            count[block[set][way].cpu]++;
        }
    }

    if (count[0] >= partitioning[pidx][0]) {    //checks if the number of block of core0 is gte to its partition
        vic_cpu = 0;    // the block that will be evicted will be from core0
    }
    else {
        vic_cpu = 1;    // the block that will be evicted will be from core1
    }

    // fill invalid line first
    for (way = 0; way < NUM_WAY; way++) {
        if (count[vic_cpu] < partitioning[pidx][vic_cpu]) {
            if ((block[set][way].valid == false)) {//If the block is invalid

                DP(if (warmup_complete[vic_cpu]) {
                    cout << "[" << NAME << "] " << func << " instr_id: " << instr_id << " invalid set: " << set << " way: " << way;
                    cout << hex << " address: " << (full_addr >> LOG2_BLOCK_SIZE) << " victim address: " << block[set][way].address << " data: " << block[set][way].data;
                    cout << dec << " lru: " << block[set][way].lru << endl;
                });

                break;
            }
        }
        else {
            break;
        }
    }

    // LRU victim

    if (way == NUM_WAY) {
        for (way = 0; way < NUM_WAY; way++) {
            if ((block[set][way].lru == partitioning[pidx][vic_cpu] - 1) && (block[set][way].cpu == vic_cpu)) {//evicting extra blocks belonging to a core
                DP(if (warmup_complete[cpu]) {
                    cout << "[" << NAME << "] " << func << " instr_id: " << instr_id << " replace set: " << set << " way: " << way;
                    cout << hex << " address: " << (full_addr >> LOG2_BLOCK_SIZE) << " victim address: " << block[set][way].address << " data: " << block[set][way].data;
                    cout << dec << " lru: " << block[set][way].lru << endl;
                });

                break;
            }
        }
    }
    if (way == NUM_WAY) {
        cerr << "[" << NAME << "] " << " no victim! set: " << set << endl;
        assert(0);
    }
    return way;
}

void CACHE::lru_update(uint32_t set, uint32_t way)
{
    // update lru replacement state
    for (uint32_t i = 0; i < NUM_WAY; i++) {
        if (block[set][i].lru < block[set][way].lru) {
            block[set][i].lru++;
        }
    }
    block[set][way].lru = 0; // promote to the MRU position
}

void CACHE::llc_lru_update(uint32_t set, uint32_t way, uint32_t cpu)
{
    // update lru replacement state
    for (uint32_t i = 0; i < NUM_WAY; i++) {
        if ((block[set][i].cpu == vic_cpu) && (block[set][i].lru < block[set][way].lru)) {
            block[set][i].lru++;
        }
    }
    block[set][way].lru = 0; // promote to the MRU position
}

void CACHE::replacement_final_stats()
{

}

#ifdef NO_CRC2_COMPILE
void InitReplacementState()
{

}

uint32_t GetVictimInSet(uint32_t cpu, uint32_t set, const BLOCK* current_set, uint64_t PC, uint64_t paddr, uint32_t type)
{
    return 0;
}

void UpdateReplacementState(uint32_t cpu, uint32_t set, uint32_t way, uint64_t paddr, uint64_t PC, uint64_t victim_addr, uint32_t type, uint8_t hit)
{

}

void PrintStats_Heartbeat()
{

}

void PrintStats()
{

}
#endif
