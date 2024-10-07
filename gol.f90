program gameoflife
use mpi

! Variables
integer, parameter :: nrows = 20, ncols = 20
integer :: N = 120
logical :: game(nrows + 2, ncols + 2), next_game(nrows + 2, ncols + 2)
integer :: i, j, step, mini_i, mini_j, cells_alive, true_count
integer :: ts, tf, rate
real :: rand, total_time
integer :: ierr, rank, size

! Setup MPI
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

! Initialize grid to false
game = .FALSE.

! Initialize grid points
if (rank .eq. 0) then
    game(2, 4) = .TRUE.
    game(3, 2) = .TRUE.
    game(3, 4) = .TRUE.
    game(4, 3) = .TRUE.
    game(4, 4) = .TRUE.

    ! Random cell initialization
    !do i = 1, nrows + 2
    !    do j = 1, ncols + 2
    !        call random_number(rand)

    !        if (rand < 0.5) then
    !            game(i, j) = .FALSE.
    !        else
    !            game(i, j) = .TRUE.
    !        end if
    !    end do
    !end do
endif

! Set the rows
game(1, :) = game(nrows + 1, :)
game(nrows + 2, :) = game(2, :)

! Set the columns
game(:, 1) = game(:, ncols + 1)
game(:, ncols + 2) = game(:, 2)

! Set the corners
game(1, 1) = game(nrows + 1, ncols + 1)
game(1, ncols + 2) = game(ncols + 1, 2)
game(nrows + 2, 1) = game(2, ncols + 1)
game(nrows + 2, ncols + 2) = game(2, 2)

! Count the initial number of alive cells
cells_alive = count(game(2:nrows+1, 2:ncols+1))

! Print out initial matrix and contents
if (rank .eq. 0) then
    print *, 'Initial Boolean Matrix'
    print *, 'Initial No. of Cells Alive: ', cells_alive

    do i = 2, nrows + 1
        do j = 2, ncols + 1
            if (game(i, j)) then
                write(*, '(A)', ADVANCE='NO') ' T'
            else
                write(*, '(A)', ADVANCE='NO') ' F'
            end if
        
            if (j == ncols + 1) then
                write(*,*)
            end if
        end do
    end do
end if

! Print blank line, set the number of steps, and matrix copy
if (rank .eq. 0) then
    print *, ''
    next_game = game
end if

! Performs N steps with a timer
call system_clock(count_rate=rate) 
call system_clock(ts)
do step = 1, N
    ! Broadcasts the whole grid to the other processes
    call MPI_BCAST(game, (nrows + 2) * (ncols + 2), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    do i = 2, nrows + 1
        do j = 2, ncols + 1
            true_count = 0    
            
            ! Checks all the cells around the current one
            do mini_i = i - 1, i + 1
                do mini_j = j - 1, j + 1
                    if (mini_i /= i .or. mini_j /= j) then
                        if (game(mini_i, mini_j)) then
                            true_count = true_count + 1
                        end if
                    end if
                end do
            end do
            
            ! Game of Life Rules
            if (true_count .eq. 3) then
                next_game(i, j) = .TRUE.
            else if (true_count .eq. 2) then
                continue
            else
                next_game(i, j) = .FALSE.
            end if
        end do
    end do
    
    ! Set the rows
    next_game(1, :) = next_game(nrows + 1, :)
    next_game(nrows + 2, :) = next_game(2, :)

    ! Set the columns
    next_game(:, 1) = next_game(:, ncols + 1)
    next_game(:, ncols + 2) = next_game(:, 2)

    ! Set the corners
    next_game(1, 1) = next_game(nrows + 1, ncols + 1)
    next_game(1, ncols + 2) = next_game(ncols + 1, 2)
    next_game(nrows + 2, 1) = next_game(2, ncols + 1)
    next_game(nrows + 2, ncols + 2) = next_game(2, 2)

    ! Gathers all the processes together
    call MPI_GATHER(next_game, (nrows + 2) * (ncols + 2) / size, &
                    MPI_INTEGER, grid, (nrows + 2)*(ncols + 2)/size, &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) then
        game = next_game
    end if
end do
call system_clock(tf)
total_time = tf - ts

! Count the final number of alive cells
cells_alive = count(game(2:nrows+1, 2:ncols+1))

! Prints out resulting grid
if (rank .eq. 0) then
    print *, 'Grid After ', N, ' steps'
    print *, 'Total Time: ', total_time
    print *, 'Final No. of Cells Alive: ', cells_alive


    do i = 2, nrows + 1
        do j = 2, ncols + 1
            if (next_game(i, j)) then
                write(*, '(A)', ADVANCE='NO') ' T'
            else
                write(*, '(A)', ADVANCE='NO') ' F'
            end if
        
            if (j == ncols + 1) then
                write(*,*)
            end if
        end do
    end do
end if

call MPI_FINALIZE(ierr)
end program
