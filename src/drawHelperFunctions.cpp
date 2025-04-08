#include <drawHelperFunctions.h>

void updateV(Eigen::MatrixXd& V, Eigen::VectorXd& q){
	int n = V.rows();  

	for (int i = 0; i < n; i++){
		V.row(i) = q.segment<3>(3 * i);
	}
    
}