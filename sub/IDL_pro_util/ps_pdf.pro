pro PS_PDF, REMOVE=REMOVE

spawn,'ls -1 *.ps > conv'
BATCH_CONV,'conv'
spawn,'source conv_'
if keyword_set(REMOVE) then spawn,'rm *.ps'
spawn,'rm conv'
spawn,'rm conv_'

return
end
