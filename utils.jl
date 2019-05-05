function string_to_array(s::String)
    temp_array = Array{Char}(UndefInitializer(),length(s))
    for i in 1:1:length(s)
        temp_array[i] = s[i]
    end
    return parse.(Type,temp_array)
end


