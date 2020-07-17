#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#ifdef MPI
#include <mpi.h>
#endif

#ifdef OPENMP
#include <omp.h>
#endif

#define ALIVE 1
#define DEAD 0

const int MINIMUM_ROWS = 1;
const int MINIMUM_COLUMNS = 1;
const int MINIMUM_TIME_STEPS = 1;
void pluralize_value_if_needed(int value)
{
    if(value != 1)
        fprintf(stderr, "s");
    return;
}
int assert_minimum_value(char which_value[16], int actual_value,
        int expected_value)
{
    int retval;

    if(actual_value < expected_value)
    {
        fprintf(stderr, "ERROR: %d %s", actual_value, which_value);
        pluralize_value_if_needed(actual_value);
        fprintf(stderr, "; need at least %d %s", expected_value, which_value);
        pluralize_value_if_needed(expected_value);
        fprintf(stderr, "\n");
        retval = -1;
    }
    else
        retval = 0;

    return retval;
}
void exit_if(int boolean_expression, char function_name[32], int OUR_RANK)
{
    if(boolean_expression)
    {
#ifdef MPI
        fprintf(stderr, "Rank %d ", OUR_RANK);
#endif
#ifdef OPENMP
        fprintf(stderr, "Thread %d ", omp_get_thread_num());
#endif
        fprintf(stderr, "ERROR in %s\n", function_name);
        exit(-1);
    }

    return;
}
int main(int argc, char **argv)
{
	
    int NUMBER_OF_ROWS = 5, NUMBER_OF_COLUMNS = 5, NUMBER_OF_TIME_STEPS = 5, 
        OUR_NUMBER_OF_ROWS = 5, OUR_RANK = 0, NUMBER_OF_PROCESSES = 1, 
        our_current_row, my_current_column, my_neighbor_row, my_neighbor_column,
        my_number_of_alive_neighbors, c, return_value, next_lowest_rank, 
        next_highest_rank;

    int **our_current_grid, **our_next_grid;
    int current_time_step;
#ifdef SHOW_RESULTS
    int current_rank;

#endif
#ifdef MPI
    double start = MPI_Wtime();
    exit_if((MPI_Init(&argc, &argv) != MPI_SUCCESS), "MPI_Init", OUR_RANK);
    exit_if((MPI_Comm_rank(MPI_COMM_WORLD, &OUR_RANK) != MPI_SUCCESS),
            "MPI_Comm_rank", OUR_RANK);
    exit_if((MPI_Comm_size(MPI_COMM_WORLD, &NUMBER_OF_PROCESSES)
                != MPI_SUCCESS), "MPI_Comm_size", OUR_RANK);
#endif
    while((c = getopt(argc, argv, "r:c:t:")) != -1)
    {
        switch(c)
        {
            case 'r':
                NUMBER_OF_ROWS = atoi(optarg);
                break;
            case 'c':
                NUMBER_OF_COLUMNS = atoi(optarg);
                break;
            case 't':
                NUMBER_OF_TIME_STEPS = atoi(optarg);
                break;
            case '?':
            default:
#ifdef MPI
                fprintf(stderr, "Usage: mpirun -np NUMBER_OF_PROCESSES %s [-r NUMBER_OF_ROWS] [-c NUMBER_OF_COLUMNS] [-t NUMBER_OF_TIME_STEPS]\n", argv[0]);
#else
                fprintf(stderr, "Usage: %s [-r NUMBER_OF_ROWS] [-c NUMBER_OF_COLUMNS] [-t NUMBER_OF_TIME_STEPS]\n", argv[0]);
#endif
                exit(-1);
        }
    }
    argc -= optind;
    argv += optind;
    return_value = assert_minimum_value("row", NUMBER_OF_ROWS, MINIMUM_ROWS);
    return_value += assert_minimum_value("column", NUMBER_OF_COLUMNS,
            MINIMUM_COLUMNS);
    return_value += assert_minimum_value("time step", NUMBER_OF_TIME_STEPS,
            MINIMUM_TIME_STEPS);
    if(return_value != 0)
        exit(-1);
    OUR_NUMBER_OF_ROWS = NUMBER_OF_ROWS / NUMBER_OF_PROCESSES;
    if(OUR_RANK == NUMBER_OF_PROCESSES - 1)
    {
        OUR_NUMBER_OF_ROWS += NUMBER_OF_ROWS % NUMBER_OF_PROCESSES;
    }
    exit_if(((our_current_grid = (int**)malloc((OUR_NUMBER_OF_ROWS + 2) 
                        * (NUMBER_OF_COLUMNS + 2) * sizeof(int))) == NULL),
            "malloc(our_current_grid)", OUR_RANK);
    exit_if(((our_next_grid = (int**)malloc((OUR_NUMBER_OF_ROWS + 2) 
                        * (NUMBER_OF_COLUMNS + 2) * sizeof(int))) == NULL),
            "malloc(our_next_grid)", OUR_RANK);
    for(our_current_row = 0; our_current_row <= OUR_NUMBER_OF_ROWS + 1;
            our_current_row++)
    {
        exit_if(((our_current_grid[our_current_row]
                        = (int*)malloc((NUMBER_OF_COLUMNS + 2) * sizeof(int))) 
                    == NULL), "malloc(our_current_grid[some_row])", OUR_RANK);
        exit_if(((our_next_grid[our_current_row]
                        = (int*)malloc((NUMBER_OF_COLUMNS + 2) * sizeof(int)))
                    == NULL), "malloc(our_next_grid[some_row])", OUR_RANK);
    }
    for(our_current_row = 1; our_current_row <= OUR_NUMBER_OF_ROWS;
            our_current_row++)
    {
#pragma omp parallel for private(my_current_column)
        for(my_current_column = 1; my_current_column <= NUMBER_OF_COLUMNS;
                my_current_column++)
        {
            our_current_grid[our_current_row][my_current_column] =
                random() % (ALIVE + 1);
        }
    }
    if(OUR_RANK == 0)
        next_lowest_rank = NUMBER_OF_PROCESSES - 1;
    else
        next_lowest_rank = OUR_RANK - 1;
    if(OUR_RANK == NUMBER_OF_PROCESSES - 1)
        next_highest_rank = 0;
    else
        next_highest_rank = OUR_RANK + 1;

    for(current_time_step = 0; current_time_step <= NUMBER_OF_TIME_STEPS - 1;
            current_time_step++)
    {
#ifdef MPI
        exit_if((MPI_Send(our_current_grid[1], NUMBER_OF_COLUMNS + 2,
                        MPI_INT, next_lowest_rank, 0, MPI_COMM_WORLD) !=
                    MPI_SUCCESS),
                "MPI_Send(top row)", OUR_RANK);

        exit_if((MPI_Send(our_current_grid[OUR_NUMBER_OF_ROWS], 
                        NUMBER_OF_COLUMNS + 2, MPI_INT, next_highest_rank,
                        0, MPI_COMM_WORLD) != MPI_SUCCESS),
                "MPI_Send(bottom row)", OUR_RANK);

        exit_if((MPI_Recv(our_current_grid[OUR_NUMBER_OF_ROWS + 1],
                        NUMBER_OF_COLUMNS + 2, MPI_INT, next_highest_rank, 
                        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) 
                    != MPI_SUCCESS),
                "MPI_Recv(bottom row)", OUR_RANK);

        exit_if((MPI_Recv(our_current_grid[0], NUMBER_OF_COLUMNS + 2,
                        MPI_INT, next_lowest_rank, 0, MPI_COMM_WORLD,
                        MPI_STATUS_IGNORE) != MPI_SUCCESS),
                "MPI_Recv(top row)", OUR_RANK);
#else
#pragma omp parallel private(my_current_column)
        for(my_current_column = 0;
                my_current_column <= NUMBER_OF_COLUMNS + 1;
                my_current_column++)
        {
            our_current_grid[0][my_current_column]
                = our_current_grid[OUR_NUMBER_OF_ROWS][my_current_column];

            our_current_grid[OUR_NUMBER_OF_ROWS + 1][my_current_column]
                = our_current_grid[1][my_current_column];
        }
#endif

        for(our_current_row = 0; our_current_row <= OUR_NUMBER_OF_ROWS + 1;
                our_current_row++)
        {
            our_current_grid[our_current_row][0] =
                our_current_grid[our_current_row][NUMBER_OF_COLUMNS];

            our_current_grid[our_current_row][NUMBER_OF_COLUMNS + 1] =
                our_current_grid[our_current_row][1];
        }

#ifdef SHOW_RESULTS
        for(current_rank = 0; current_rank <= NUMBER_OF_PROCESSES - 1;
                current_rank++)
            printf("Time Step %d, Rank %d:\n", current_time_step, OUR_RANK);
        for(our_current_row = 0; our_current_row <= OUR_NUMBER_OF_ROWS + 1;
                our_current_row++)
        {
            if(our_current_row == 1)
            {
                for(my_current_column = 0;
                        my_current_column <= NUMBER_OF_COLUMNS + 1 + 2;
                        my_current_column++)
                {
                    printf("- ");
                }
                printf("\n");
            }

            for(my_current_column = 0;
                    my_current_column <= NUMBER_OF_COLUMNS + 1; 
                    my_current_column++)
            {
                if(my_current_column == 1)
                {
                    printf("| ");
                }

                printf("%d ", our_current_grid[our_current_row]
                        [my_current_column]);

                if(my_current_column == NUMBER_OF_COLUMNS)
                {
                    printf("| ");
                }
            }
            printf("\n");

            if(our_current_row == OUR_NUMBER_OF_ROWS)
            {
                for(my_current_column = 0;
                        my_current_column <= NUMBER_OF_COLUMNS + 1 + 2;
                        my_current_column++)
                {
                    printf("- ");
                }
                printf("\n");
            }
        }
#endif

        for(our_current_row = 1; our_current_row <= OUR_NUMBER_OF_ROWS;
                our_current_row++)
        {
#pragma omp parallel for private(my_current_column, my_neighbor_row, my_neighbor_column, my_number_of_alive_neighbors)
            for(my_current_column = 1; my_current_column <= NUMBER_OF_COLUMNS;
                    my_current_column++)
            {
                my_number_of_alive_neighbors = 0;

                for(my_neighbor_row = our_current_row - 1;
                        my_neighbor_row <= our_current_row + 1;
                        my_neighbor_row++)
                {

                    for(my_neighbor_column = my_current_column - 1;
                            my_neighbor_column <= my_current_column + 1;
                            my_neighbor_column++)
                    {

                        if((my_neighbor_row != our_current_row
                                    || my_neighbor_column != my_current_column)
                                && (our_current_grid[my_neighbor_row]
                                    [my_neighbor_column] == ALIVE))
                        {

                            my_number_of_alive_neighbors++;
                        }
                    }
                }

                if(my_number_of_alive_neighbors < 2)
                {
                    our_next_grid[our_current_row][my_current_column] = DEAD;
                }

                if(our_current_grid[our_current_row][my_current_column] == ALIVE
                        && (my_number_of_alive_neighbors == 2
                            || my_number_of_alive_neighbors == 3))
                {
                    our_next_grid[our_current_row][my_current_column] = ALIVE;
                }

                if(my_number_of_alive_neighbors > 3)
                {
                    our_next_grid[our_current_row][my_current_column] = DEAD;
                }


                if(our_current_grid[our_current_row][my_current_column] == DEAD
                        && my_number_of_alive_neighbors == 3)
                {
                    our_next_grid[our_current_row][my_current_column] = ALIVE;
                }
            }
        }


        for(our_current_row = 1; our_current_row <= OUR_NUMBER_OF_ROWS;
                our_current_row++)
        {
#pragma omp parallel for private(my_current_column)
            for(my_current_column = 1; my_current_column <= NUMBER_OF_COLUMNS;
                    my_current_column++)
            {
                our_current_grid[our_current_row][my_current_column] =
                    our_next_grid[our_current_row][my_current_column];
            }
        }
    }


    for(our_current_row = OUR_NUMBER_OF_ROWS + 1; our_current_row >= 0;
            our_current_row--)
    {
        free(our_next_grid[our_current_row]);
        free(our_current_grid[our_current_row]);
    }
    free(our_next_grid);
    free(our_current_grid);


#ifdef MPI
	double end = MPI_Wtime();
	printf("Total time is: %f",end-start);
    exit_if((MPI_Finalize() != MPI_SUCCESS), "MPI_Finalize", OUR_RANK);
#endif
	
    return 0;
}
