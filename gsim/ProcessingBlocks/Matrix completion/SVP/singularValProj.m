function CompletedMatSVP =singularValProj(X,o,desiredRank)
mc = MatrixCompletor_SVP;
mc.alpha_init= 1;
mc.desiredRank= desiredRank;
mc.gamma_const = 5e-1;
mc.mask = o;
mc.X_0 =X;
mc.delta_const = 0.8;
[X_out,~]=mc.singularValueProjection();
CompletedMatSVP=X_out;
end

