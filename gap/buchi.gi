xxx_global := 42;

Print("Loading buchi functionality...\n");

InstallMethod(BuchiMachine, "(FR) for a Julia Büchi machine", [IsJuliaObject],
              function(obj)
                  local F, M;
                  F := FRMFamily([1..10]);
                  M := Objectify(NewType(F, IsMealyElement and IsBuchiMachineRep and IsJuliaWrapper), rec());
                  SetJuliaPointer(M,obj);
                  return M;
              end);
