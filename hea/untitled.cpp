double * alpha, * beta;
		alpha = (double *) calloc(M+1, sizeof(double));
		beta = (double *) calloc(M+1, sizeof(double));

		for ( int t = 0; t < P; t++ ){
			tmp = prevLayer;
			prevLayer = currLayer;
			currLayer = tmp;
			int x;
			alpha[0] = 0;
			beta[0] = 0;
			int xi1 = a1(t)/prevLayer[1];
			int mu1 = a1(t) - prevLayer[1]*xi1;
			alpha[1] = xi1;
			beta[1] = mu1;
			for ( x = 1; x < M; x++ ){
				alpha[x+1] = dt*a*a/(dx*dx) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
				beta[x+1] = ( dt*f(x*dx,t*dt,a) + prevLayer[x] + beta[x]*(dt*a*a)/(dx*dx) ) / ( 1 + 2*dt*a*a/(dx*dx) - alpha[x]*dt*a*a/(dx*dx) );
			}
			int xi2 = a2(t) / prevLayer[x-1];
			int mu2 = a2(t) - prevLayer[x-1]*xi2;
			currLayer[M] = (xi2*beta[x]+mu2)/(1-xi2*alpha[x]);
			for ( int x = M; x > 0/*1*/; x-- ){
				currLayer[x-1] = alpha[x]*currLayer[x]+beta[x];
			}
			
		}