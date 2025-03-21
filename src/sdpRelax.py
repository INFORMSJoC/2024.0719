

one = zeros(1, m+1);
begin_index = m;
iter_index = 1;
for i = 1:m
    for j = i+1:m
        iter_index = iter_index + 1;
    end
end

begin_index = begin_index+iter_index-1;
iter_index = 1;
for i = 1:m
    for j = i:m
        iter_index = iter_index + 1;
    end
end


begin_index = begin_index+iter_index-1;

%temp_b{begin_index+2} = -(v-g)' * mu;
one(1) = begin_index+2;

begin_index = begin_index + 2;


begin_index = begin_index + m;
for i = 1:m
    temp = -0.5*A1Complete'*(v-g);
    %temp_b{begin_index+i} = temp(i);
    one(1+i) = begin_index+i;
end

two = ones(1, m+1);

three = zeros(1, m+1);
three(1) = -(v-g)' * mu;
for i = 1:m
    temp = -0.5*A1Complete'*(v-g);
    %temp_b{begin_index+i} = temp(i);
    three(1+i) = temp(i);
end

b_input = sparse(one, two, three, numberOfConstraints, 1);



% A_input


% for i = 1:m
%     temp_A{i}(K_input.l + 1) = 1;
%     temp_A{i}(K_input.l + (m+1)*i + 1+i) = -1;
% end

nonzeroNumber = 2*m;

begin_index = m;
iter_index = 1;
for i = 1:m
    for j = i+1:m
        % temp_A{begin_index+iter_index}(K_input.l + (m+1)*i + 1+j) = 1;
        nonzeroNumber = nonzeroNumber + 1;
        iter_index = iter_index + 1;
    end
end

begin_index = begin_index+iter_index-1;
iter_index = 1;
for i = 1:m
    for j = i:m
%         temp_A{begin_index+iter_index}(K_input.l + (m+1)^2 + (m+1)*i + 1+j) = 1;
%         temp_A{begin_index+iter_index}(K_input.l + 2*(m+1)^2 + (m+1)*i + 1+j) = -1;
        nonzeroNumber = nonzeroNumber + 2;
        iter_index = iter_index + 1;
    end
end


begin_index = begin_index+iter_index-1;


% temp_A{begin_index+1}(K_input.l+ (m+1)^2 + 1) = -1;
% temp_A{begin_index+1}(K_input.l-1) = 1;
% temp_A{begin_index+1}(K_input.l) = -1;
% for i = 1:n
%     temp_A{begin_index+1}(i) = -(c(i)-v(i)); 
% end
% for i = 1:2*m
%     temp = A*mu-b;
%     temp_A{begin_index+1}(n+i) = temp(i);
% end
nonzeroNumber = nonzeroNumber + 3 + n + 2*m;



% temp_A{begin_index+2}(K_input.l+ 2*(m+1)^2 + 1) = -1;
% temp_A{begin_index+2}(K_input.l-1) = 1;
% temp_A{begin_index+2}(K_input.l) = -1;
% for i = 1:n
%     temp_A{begin_index+2}(i) = -(c(i)-g(i)); 
% end
% for i = 1:2*m
%     temp = A*mu-b;
%     temp_A{begin_index+2}(n+2*m+i) = temp(i);
% end
nonzeroNumber = nonzeroNumber + 3 + n + 2*m;

begin_index = begin_index + 2;


for i = 1:m
%     temp_A{begin_index+i}(K_input.l + (m+1)^2 + 1 + i) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
    nonzeroNumber = nonzeroNumber + 2;
    temp = 0.5*A1Complete'*A';
    temp = temp(i,:);
    for j = 1:2*m
%         temp_A{begin_index+i}(n+j) = temp(j);
        nonzeroNumber = nonzeroNumber + 1;
    end
end

begin_index = begin_index + m;
for i = 1:m
%     temp_A{begin_index+i}(K_input.l + 2*(m+1)^2 + 1 + i) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
    nonzeroNumber = nonzeroNumber + 2;
    temp = 0.5*A1Complete'*A';
    temp = temp(i,:);
    for j = 1:2*m
%         temp_A{begin_index+i}(n+2*m+j) = temp(j);
        nonzeroNumber = nonzeroNumber + 1;
    end
end


one = zeros(1, nonzeroNumber);
two = zeros(1, nonzeroNumber);
three = zeros(1, nonzeroNumber);


for i = 1:m
%     temp_A{i}(K_input.l + 1) = 1;
%     temp_A{i}(K_input.l + (m+1)*i + 1+i) = -1;
    one(2*i-1) = i;
    two(2*i-1) = K_input.l + 1;
    three(2*i-1) = 1;

    one(2*i) = i;
    two(2*i) = K_input.l + (m+1)*i + 1+i;
    three(2*i) = -1;
end

nonzeroNumber = 2*m;


begin_index = m;
iter_index = 1;
for i = 1:m
    for j = i+1:m
        % temp_A{begin_index+iter_index}(K_input.l + (m+1)*i + 1+j) = 1;
        nonzeroNumber = nonzeroNumber + 1;
        one(nonzeroNumber) = begin_index+iter_index;
        two(nonzeroNumber) = K_input.l + (m+1)*i + 1+j;
        three(nonzeroNumber) = 1;
        iter_index = iter_index + 1;
    end
end

begin_index = begin_index+iter_index-1;
iter_index = 1;
for i = 1:m
    for j = i:m
%         temp_A{begin_index+iter_index}(K_input.l + (m+1)^2 + (m+1)*i + 1+j) = 1;
%         temp_A{begin_index+iter_index}(K_input.l + 2*(m+1)^2 + (m+1)*i + 1+j) = -1;
        nonzeroNumber = nonzeroNumber + 2;

        one(nonzeroNumber-1) = begin_index+iter_index;
        two(nonzeroNumber-1) = K_input.l + (m+1)^2 + (m+1)*i + 1+j;
        three(nonzeroNumber-1) = 1;

        one(nonzeroNumber) = begin_index+iter_index;
        two(nonzeroNumber) = K_input.l + 2*(m+1)^2 + (m+1)*i + 1+j;
        three(nonzeroNumber) = -1;
        iter_index = iter_index + 1;
    end
end


begin_index = begin_index+iter_index-1;


% temp_A{begin_index+1}(K_input.l+ (m+1)^2 + 1) = -1;
one(nonzeroNumber+1) = begin_index+1;
two(nonzeroNumber+1) = K_input.l+ (m+1)^2 + 1;
three(nonzeroNumber+1) = -1;
% temp_A{begin_index+1}(K_input.l-1) = 1;
one(nonzeroNumber+2) = begin_index+1;
two(nonzeroNumber+2) = K_input.l-1;
three(nonzeroNumber+2) = 1;
% temp_A{begin_index+1}(K_input.l) = -1;
one(nonzeroNumber+3) = begin_index+1;
two(nonzeroNumber+3) = K_input.l;
three(nonzeroNumber+3) = -1;
for i = 1:n
%     temp_A{begin_index+1}(i) = -(c(i)-v(i)); 
    one(nonzeroNumber+3+i) = begin_index+1;
    two(nonzeroNumber+3+i) = i;
    three(nonzeroNumber+3+i) = -(c(i)-v(i));
end
for i = 1:2*m
    temp = A*mu-b;
%     temp_A{begin_index+1}(n+i) = temp(i);
    one(nonzeroNumber+3+n+i) = begin_index+1;
    two(nonzeroNumber+3+n+i) = n+i;
    three(nonzeroNumber+3+n+i) = temp(i);
end
nonzeroNumber = nonzeroNumber + 3 + n + 2*m;



% temp_A{begin_index+2}(K_input.l+ 2*(m+1)^2 + 1) = -1;
one(nonzeroNumber+1) = begin_index+2;
two(nonzeroNumber+1) = K_input.l+ 2*(m+1)^2 + 1;
three(nonzeroNumber+1) = -1;
% temp_A{begin_index+2}(K_input.l-1) = 1;
one(nonzeroNumber+2) = begin_index+2;
two(nonzeroNumber+2) = K_input.l-1;
three(nonzeroNumber+2) = 1;
% temp_A{begin_index+2}(K_input.l) = -1;
one(nonzeroNumber+3) = begin_index+2;
two(nonzeroNumber+3) = K_input.l;
three(nonzeroNumber+3) = -1;
for i = 1:n
%     temp_A{begin_index+2}(i) = -(c(i)-g(i)); 
    one(nonzeroNumber+3+i) = begin_index+2;
    two(nonzeroNumber+3+i) = i;
    three(nonzeroNumber+3+i) = -(c(i)-g(i));
end
for i = 1:2*m
    temp = A*mu-b;
%     temp_A{begin_index+2}(n+2*m+i) = temp(i);
    one(nonzeroNumber+3+n+i) = begin_index+2;
    two(nonzeroNumber+3+n+i) = n+2*m+i;
    three(nonzeroNumber+3+n+i) = temp(i);
end
nonzeroNumber = nonzeroNumber + 3 + n + 2*m;

begin_index = begin_index + 2;


for i = 1:m
%     temp_A{begin_index+i}(K_input.l + (m+1)^2 + 1 + i) = -1;
    one(nonzeroNumber+1) = begin_index+i;
    two(nonzeroNumber+1) = K_input.l + (m+1)^2 + 1 + i;
    three(nonzeroNumber+1) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
    one(nonzeroNumber+2) = begin_index+i;
    two(nonzeroNumber+2) = K_input.l + 1 + i;
    three(nonzeroNumber+2) = 0.5;
    nonzeroNumber = nonzeroNumber + 2;
    temp = 0.5*A1Complete'*A';
    temp = temp(i,:);
    for j = 1:2*m
%         temp_A{begin_index+i}(n+j) = temp(j);
        nonzeroNumber = nonzeroNumber + 1;
        one(nonzeroNumber) = begin_index+i;
        two(nonzeroNumber) = n+j;
        three(nonzeroNumber) = temp(j);
    end
end

begin_index = begin_index + m;
for i = 1:m
%     temp_A{begin_index+i}(K_input.l + 2*(m+1)^2 + 1 + i) = -1;
    one(nonzeroNumber+1) = begin_index+i;
    two(nonzeroNumber+1) = K_input.l + 2*(m+1)^2 + 1 + i;
    three(nonzeroNumber+1) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
    one(nonzeroNumber+2) = begin_index+i;
    two(nonzeroNumber+2) = K_input.l + 1 + i;
    three(nonzeroNumber+2) = 0.5;
    nonzeroNumber = nonzeroNumber + 2;
    temp = 0.5*A1Complete'*A';
    temp = temp(i,:);
    for j = 1:2*m
%         temp_A{begin_index+i}(n+2*m+j) = temp(j);
        nonzeroNumber = nonzeroNumber + 1;
        one(nonzeroNumber) = begin_index+i;
        two(nonzeroNumber) = n+2*m+j;
        three(nonzeroNumber) = temp(j);
    end
end


A_input = sparse(one, two, three, numberOfConstraints, X_dim);



% whos("c_input", "b_input", "A_input");
% return;

% temp_A = {};
% temp_b = {};
% 
% 
% for i = 1:m
%     temp_A{i} = sparse(1, X_dim);
%     temp_A{i}(K_input.l + 1) = 1;
%     temp_A{i}(K_input.l + (m+1)*i + 1+i) = -1;
%     temp_b{i} = 0;
% end
% 
% begin_index = m;
% iter_index = 1;
% for i = 1:m
%     for j = i+1:m
%         temp_A{begin_index+iter_index} = sparse(1,X_dim);
%         temp_A{begin_index+iter_index}(K_input.l + (m+1)*i + 1+j) = 1;
%         temp_b{begin_index+iter_index} = 0;
%         iter_index = iter_index + 1;
%     end
% end
% 
% begin_index = begin_index+iter_index-1;
% iter_index = 1;
% for i = 1:m
%     for j = i:m
%         temp_A{begin_index+iter_index} = sparse(1,X_dim);
%         temp_A{begin_index+iter_index}(K_input.l + (m+1)^2 + (m+1)*i + 1+j) = 1;
%         temp_A{begin_index+iter_index}(K_input.l + 2*(m+1)^2 + (m+1)*i + 1+j) = -1;
%         temp_b{begin_index+iter_index} = 0;
%         iter_index = iter_index + 1;
%     end
% end
% 
% 
% begin_index = begin_index+iter_index-1;
% 
% temp_A{begin_index+1} = sparse(1,X_dim);
% temp_A{begin_index+1}(K_input.l+ (m+1)^2 + 1) = -1;
% temp_A{begin_index+1}(K_input.l-1) = 1;
% temp_A{begin_index+1}(K_input.l) = -1;
% for i = 1:n
%     temp_A{begin_index+1}(i) = -(c(i)-v(i)); 
% end
% for i = 1:2*m
%     temp = A*mu-b;
%     temp_A{begin_index+1}(n+i) = temp(i);
% end
% temp_b{begin_index+1} = 0;
% 
% 
% temp_A{begin_index+2} = sparse(1,X_dim);
% temp_A{begin_index+2}(K_input.l+ 2*(m+1)^2 + 1) = -1;
% temp_A{begin_index+2}(K_input.l-1) = 1;
% temp_A{begin_index+2}(K_input.l) = -1;
% for i = 1:n
%     temp_A{begin_index+2}(i) = -(c(i)-g(i)); 
% end
% for i = 1:2*m
%     temp = A*mu-b;
%     temp_A{begin_index+2}(n+2*m+i) = temp(i);
% end
% temp_b{begin_index+2} = -(v-g)' * mu;
% 
% begin_index = begin_index + 2;
% 
% 
% for i = 1:m
%     temp_A{begin_index+i} = sparse(1,X_dim);
%     temp_A{begin_index+i}(K_input.l + (m+1)^2 + 1 + i) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
%     temp = 0.5*A1Complete'*A';
%     temp = temp(i,:);
%     for j = 1:2*m
%         temp_A{begin_index+i}(n+j) = temp(j);
%     end
%     temp_b{begin_index+i} = 0;
% end
% 
% begin_index = begin_index + m;
% for i = 1:m
%     temp_A{begin_index+i} = sparse(1,X_dim);
%     temp_A{begin_index+i}(K_input.l + 2*(m+1)^2 + 1 + i) = -1;
%     temp_A{begin_index+i}(K_input.l + 1 + i) = 0.5;
%     temp = 0.5*A1Complete'*A';
%     temp = temp(i,:);
%     for j = 1:2*m
%         temp_A{begin_index+i}(n+2*m+j) = temp(j);
%     end
%     temp = -0.5*A1Complete'*(v-g);
%     temp_b{begin_index+i} = temp(i);
% end
% 
% numberOfConstraints = begin_index+m;
% 
% 
% 
% 
% A_input = sparse(numberOfConstraints, X_dim);
% b_input = sparse(numberOfConstraints, 1);
% for i = 1:numberOfConstraints
%     A_input(i,:) = temp_A{i};
%     b_input(i) = temp_b{i};
% end



tic
[X_opt, n_y, n_info, n_r] = sdplr(A_input, b_input, c_input, K_input);
f_opt = c_input' * X_opt;
 
CPUTime = toc;



% c_input1 = c_input(1:n+4*m+2);
% c_input2 = reshape(c_input(n+4*m+3:n+4*m+2+(m+1)^2), m+1, m+1);
% c_input3 = reshape(c_input(n+4*m+2+(m+1)^2+1:n+4*m+2+2*(m+1)^2), m+1, m+1);
% c_input4 = reshape(c_input(n+4*m+2+2*(m+1)^2+1:n+4*m+2+3*(m+1)^2), m+1, m+1);
% tic 
% 
% cvx_begin sdp
%  
%                 variable X1(n+4*m+2) nonnegative;
%                 variable X2(m+1, m+1) symmetric;
%                 variable X3(m+1, m+1) symmetric;
%                 variable X4(m+1, m+1) symmetric;
%                                 
%                 %%
%                 minimize( c_input1' * X1 + sum(sum( c_input2.*X2 )) + sum(sum( c_input3.*X3 )) + sum(sum( c_input4.*X4 )) );
%                 %%
%                 subject to
%                 
%                 for i = 1:numberOfConstraints
%                     temp_A1 = A_input(i,1:n+4*m+2);
%                     temp_A2 = reshape(A_input(i, n+4*m+3:n+4*m+2+(m+1)^2), m+1, m+1);
%                     temp_A3 = reshape(A_input(i, n+4*m+2+(m+1)^2+1:n+4*m+2+2*(m+1)^2), m+1, m+1);
%                     temp_A4 = reshape(A_input(i, n+4*m+2+2*(m+1)^2+1:n+4*m+2+3*(m+1)^2), m+1, m+1);
%                     temp_A1 * X1 + sum(sum( temp_A2.*X2 )) + sum(sum( temp_A3.*X3 )) + sum(sum( temp_A4.*X4 )) == b_input(i);
%                 end
%                 X2 >= 0;
%                 X3 >= 0;
%                 X4 >= 0;
%                 
%                 
%  cvx_end
%     
% CPUTime = toc;  
%  
%  
%     disp(['Problem is ' cvx_status])
%     if ~strfind(cvx_status,'Solved')
%       return
%     end
%     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
%     disp(['optimal value of cvx:',num2str(cvx_optval)]);
%     f_opt=cvx_optval;
%     X_opt=X1;

end


