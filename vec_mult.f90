module vec_mult
    implicit none

contains

    subroutine vec_mul(n, eigen_matr1, eigen_matr2, overlap, res_matr)

    ! subroutine to calculate the normalization rule for generalized eigenvalue problem, which means
    ! generalized eigenvalue problem A*x = (lambda)*B*x
    ! if Z is matrix of eigenvectors then
    ! Z**H*B*Z = I

    ! here we assume
    ! eigen_matr2**H*overlap*eigen_matr1 = res_matr

    ! input
    integer :: n
    ! n - matrix dimension
    complex :: eigen_matr1(n,n), eigen_matr2(n,n), overlap(n,n), res_matr(n,n)
    integer :: i,j
    complex :: matr1(n,n), matr2(n,n)

    do i=1,8
        do j=1,8
            matr2(i,j) = conjg(eigen_matr2(j,i))
        enddo
    enddo

    matr1 = matmul(overlap, eigen_matr1)
    res_matr = matmul(matr2,matr1)

    return

    end subroutine vec_mul

end module vec_mult
