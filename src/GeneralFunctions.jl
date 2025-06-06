function bold_integer(i::Int)::String 
  dgts = digits( i, base = 10 )
  last( cumprod( [ bold_digits_dict[a] for a in dgts ] ) )
end 

bold_digits_dict = 
  Dict(
    0 => "𝟎", 1 => "𝟏", 2 => "𝟐", 3 => "𝟑", 4 => "𝟒", 
    5 => "𝟓", 6 => "𝟔", 7 => "𝟕", 8 => "𝟖", 9 => "𝟗"
  )

function is_constant_array( arr ) 
  if isempty(arr)
    return true 
  end
  first = arr[1]
  return all( element === first for element in arr )
end
