function scene = fromSnapshot(snapshot,num)
%FROMSNAPSHOT Replay saved coordinates without consuming randomness.
scene = snapshot;
assert(isequal(size(scene.bsPosM),[num.K,2]),'Snapshot BS dimension mismatch.');
assert(isequal(size(scene.pUserPosM),[num.M,2]),'Snapshot P-user dimension mismatch.');
assert(isequal(size(scene.cUserPosM),[num.K,num.N,2]),'Snapshot C-user dimension mismatch.');
assert(all(isfinite(scene.bsPosM(:))) && all(isfinite(scene.pUserPosM(:))) ...
    && all(isfinite(scene.cUserPosM(:))),'Snapshot contains invalid coordinates.');
end
