function output = log_sum_exp(v_input)

% output = log(sum(exp(v_input))) with a trick for numerical stability

M = max(v_input);
output = log(sum(exp(v_input-M)))+M;

end

