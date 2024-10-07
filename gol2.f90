program gameoflife
use mpi

! Variables
integer, parameter :: nrows = 20, ncols = 20
integer :: N = 80
logical :: game(nrows + 2, ncols + 2), next_game(nrows + 2, ncols + 2)
logical, allocatable :: local_game(:,:), local_next_game(:,:)
integer :: i, j, step, mini_i, mini_j, cells_alive, true_count, total_cells_alive
real :: rand
integer :: ierr, rank, size, left, right, up, down, local_ncols
logical, allocatable :: send_left(:), recv_left(:), send_right(:), recv_right(:)
logical, allocatable :: send_up(:), recv_up(:), send_down(:), recv_down(:)

! Setup MPI
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

! Column-wise domain decomposition starts here
local_ncols = ncols / size

! Allocating memory for grids and communication buffers
! Each of the 4 processes will use its own grid
allocate(local_game(nrows + 2, local_ncols + 2))
allocate(local_next_game(nrows + 2, local_ncols + 2))
allocate(send_left(nrows + 2))
allocate(recv_left(nrows + 2))
allocate(send_right(nrows + 2))
allocate(recv_right(nrows + 2))
allocate(send_up(local_ncols + 2))
allocate(recv_up(local_ncols + 2))
allocate(send_down(local_ncols + 2))
allocate(recv_down(local_ncols + 2))

! Initialize grid to false
local_game = .FALSE.

! Initialize grid points
if (rank .eq. 0) then
    local_game(2, 4) = .TRUE.
    local_game(3, 2) = .TRUE.
    local_game(3, 4) = .TRUE.
    local_game(4, 3) = .TRUE.
    local_game(4, 4) = .TRUE.

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

! Broadcasts the whole grid to the other processes
call MPI_BCAST(local_game, (nrows + 2) * (local_ncols + 2), MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

left = mod(rank - 1 + size, size)
right = mod(rank + 1, size)
up = rank
down = rank

! Count the initial number of alive cells
local_cells_alive = count(local_game(2:nrows+1, 2:local_ncols+1))
call MPI_REDUCE(local_cells_alive, cells_alive, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

! Print out initial matrix and contents
if (rank .eq. 0) then
    print *, 'Initial Boolean Matrix'
    print *, 'Initial No. of Cells Alive: ', cells_alive
end if

do i = 2, nrows + 1
    do j = 2, local_ncols + 1
        if (local_game(i, j)) then
            write(*, '(A)', ADVANCE='NO') ' T'
        else
            write(*, '(A)', ADVANCE='NO') ' F'
        end if
        
        if (j == local_ncols + 1) then
            write(*,*)
        end if
    end do
end do

! Print blank line, set the number of steps, and matrix copy
if (rank .eq. 0) then
    print *, ''
    next_game = game
    local_next_game = local_game
end if

! Performs N steps
do step = 1, N
    ! If there is > 1 processors, exchange info about grid ghost cells
    if (size > 1) then
        send_left = local_game(:, 2)
        send_right = local_game(:, local_ncols + 1)
        send_up = local_game(2, :)
        send_down = local_game(nrows + 1, :)

        call MPI_SENDRECV(send_left, nrows + 2, MPI_LOGICAL, left, 1, &
                          recv_right, nrows + 2, MPI_LOGICAL, right, 1, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

        call MPI_SENDRECV(send_right, nrows + 2, MPI_LOGICAL, right, 2, &
                          recv_left, nrows + 2, MPI_LOGICAL, left, 2, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

        call MPI_SENDRECV(send_up, local_ncols + 2, MPI_LOGICAL, up, 3, &
                          recv_down, local_ncols + 2, MPI_LOGICAL, down, 3, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

        call MPI_SENDRECV(send_down, local_ncols + 2, MPI_LOGICAL, down, 4, &
                          recv_up, local_ncols + 2, MPI_LOGICAL, up, 4, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        local_game(:, 1) = recv_left
        local_game(:, local_ncols + 2) = recv_right
        local_game(1, :) = recv_up
        local_game(nrows + 2, :) = recv_down
    end if
    do i = 2, nrows + 1
        do j = 2, local_ncols + 1
            true_count = 0    
            
            ! Checks all the cells around the current one
            do mini_i = i - 1, i + 1
                do mini_j = j - 1, j + 1
                    if (mini_i /= i .or. mini_j /= j) then
                        if (local_game(mini_i, mini_j)) then
                            true_count = true_count + 1
                        end if
                    end if
                end do
            end do
            
            ! Game of Life Rules
            if (true_count .eq. 3) then
                local_next_game(i, j) = .TRUE.
            else if (true_count .eq. 2) then
                continue
            else
                local_next_game(i, j) = .FALSE.
            end if
        end do
    end do
    
    ! Gathers all the processes together
    call MPI_GATHER(next_game, (nrows + 2) * (ncols + 2) / size, &
                    MPI_INTEGER, grid, (nrows + 2)*(ncols + 2)/size, &
                    MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (rank .eq. 0) then
        local_game = local_next_game
        game = next_game
    end if
end do

! Count the final number of alive cells
local_cells_alive = count(local_next_game(2:nrows+1, 2:local_ncols+1))

! Prints out resulting grid
if (rank .eq. 0) then
    print *, 'Grid After ', N, ' steps'
    print *, 'Final No. of Cells Alive: ', cells_alive

    do i = 2, nrows + 1
        do j = 2, local_ncols + 1
            if (local_next_game(i, j)) then
                write(*, '(A)', ADVANCE='NO') ' T'
            else
                write(*, '(A)', ADVANCE='NO') ' F'
            end if
        
            if (j == local_ncols + 1) then
                write(*,*)
            end if
        end do
    end do
end if

call MPI_FINALIZE(ierr)
end program
