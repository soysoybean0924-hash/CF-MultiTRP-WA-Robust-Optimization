function test_target_group_interface()
% No noise/power conversion is allowed in this structural adapter.
resource=struct();
resource.H_true=complex(zeros(2,3,4,2));
resource.H_true(1,1,1,1)=2; resource.H_true(2,3,4,2)=5;
resource.servingMask=logical([1 0;1 1;0 1]);
resource.ueIds=[7;12]; resource.rb_count=16;
resource.resourceWeightNu=16/273;
Hdense=resource.H_true+1;
out=chasi_target_group_to_solver(resource,Hdense);
assert(isequal(out.serviceTRP,[1 2;2 3]));
assert(isequal(out.ueIds,[7;12]));
assert(out.rbCount==16 && out.resourceWeightNu==16/273);
assert(isequal(size(out.Htrue),[2 3]));
assert(isequal(size(out.Hdesign{1,1,1}),[4 2]));
assert(out.Htrue{2,3,1}(4,2)==5);
assert(out.Hdesign{2,3,1}(4,2)==6);
resource.servingMask(2,1)=false;
failed=false;
try
    chasi_target_group_to_solver(resource,Hdense);
catch exception
    failed=strcmp(exception.identifier,'chasi:VariableServiceCount');
end
assert(failed);
fprintf('test_target_group_interface PASS\n');
end
