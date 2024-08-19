# Manim Jr.
Some advances of formula supplements base on Janim.

## MyTypst

The `MyTypst` class extends the functionality of the `TypstDoc` class (presumably from the `typst` library) by adding methods for advanced indexing and searching capabilities. This class is particularly useful for accessing and manipulating specific parts of a document based on patterns, positions, and other criteria.

### Constructor Method: `__init__`
This method initializes the `MyTypst` object by calling the superclass constructor with the given `source` and `path`.

- **Parameters**:
  - `source`: The source of the document (e.g., file name or content).
  - `path`: The path to the source file.
  - `**kwargs`: Additional keyword arguments.

- **Behavior**:
  - Calls the superclass constructor with the provided `source` and `path`.

### Custom Indexing Method: `__getitem__`
This method overloads the indexing operator `[]` to provide flexible access to different parts of the document.

- **Parameters**:
  - `key`: The index or key to access the document.

- **Behavior**:
  - The method starts by checking the type of the key using a `match` statement.
  - If the key is an `int` or `slice`, it simply calls the superclass's `__getitem__` method to handle these basic cases.
  - If the key is a `list` containing only integers, it maps each integer to a call to the superclass's `__getitem__` method and then groups the results into a `Group` object.
  - If the key is a `list` containing only booleans, it uses the boolean values to select which elements to include in the output.
  - If the key is a `str` representing a pattern, it calls the `slices` method with the pattern and `0` as arguments to find the slice corresponding to the first occurrence of the pattern.
  - If the key is a `tuple` or `list` with a `str` followed by an `int`, it calls the `slices` method with the pattern and the ordinal number to find the slice corresponding to the nth occurrence of the pattern.
  - If the key is a `tuple` or `list` with a `str` followed by a `list` of `int`s, it calls the `slices` method with the pattern and the list of ordinal numbers to find the slices corresponding to those occurrences of the pattern.
  - If the key is a `dict`, it flattens the dictionary's values and calls the `get` method with the resulting slices.
  - If the key is a `tuple` or `list` of mixed types (patterns and ordinals), it calls the `multi_slices` method to process each pattern and ordinal combination separately and returns a `Group` of the results.
  - There is special handling for keys that represent "empty" values (like `None`, `[]`, `()`, `""`, etc.), where the entire document is returned unchanged.

### Additional Methods

#### `multi_slices`
This method finds the slices corresponding to multiple occurrences of a pattern or patterns.

- **Parameters**:
  - `patterns`: An iterable of patterns.
  - `ordinals`: An iterable of ordinal numbers indicating the nth occurrence of each pattern.

- **Behavior**:
  - It iterates over the `patterns` and `ordinals` to find the slices for each pattern and ordinal combination.

#### `slices`
This method finds the slices corresponding to the occurrence of a pattern.

- **Parameters**:
  - `pattern`: The pattern to search for.
  - `ordinal`: The ordinal number indicating the nth occurrence of the pattern.

- **Behavior**:
  - It searches the document for the specified pattern and ordinal number and returns the slice.

#### `get`
This method extracts the specified part of the document.

- **Parameters**:
  - `slice`: The slice indicating the part of the document to extract.

- **Behavior**:
  - It uses the provided slice to extract the corresponding part of the document.

### Summary
The `MyTypst` class provides a rich set of features for accessing and manipulating documents. It supports various forms of indexing, including integer indices, slices, boolean masks, and pattern-based indexing. The `__getitem__` method is the core of this functionality, allowing users to easily access specific parts of the document or search for patterns within the document. The additional methods (`multi_slices`, `slices`, `get`, and `_with_empty`) support the functionality provided by `__getitem__` and enable more complex operations on the document.

The class is designed to be user-friendly and versatile, making it suitable for a wide range of applications where precise control over document content is necessary.

## TransformInParts
This code defines a class `TransformInParts` that inherits from `AnimGroup`. The class is designed to create animations that transform items between two states, potentially involving intermediate steps like moving or transforming segments of the items. Let's break down the methods and their functionality:

### `__init__` Method
This is the constructor method for `TransformInParts`.

- **Parameters**:
  - `source`: An iterable collection of `Item` objects representing the initial state.
  - `target`: An iterable collection of `Item` objects representing the final state.
  - `durations`: An iterable of integers specifying the duration of each animation, or a single integer to be used cyclically.
  - `trs_keywords`: An iterable of dictionaries containing keyword arguments for the `Transform` animations.
  - `**kwargs`: Additional keyword arguments passed to the `AnimGroup` constructor.

- **Behavior**:
  - It processes the `trs_keywords` and `durations` to ensure they are in the correct format (iterables).
  - If `durations` is provided, it cycles through the durations and assigns them to the `trs_keywords`.
  - It creates `Transform` animations for each pair of items in `source` and `target` using the processed `trs_keywords`.
  - Finally, it initializes the `AnimGroup` with the created animations.

### `from_moved` Static Method
This method handles transformations starting from a moved state.

- **Parameters**:
  - `source`: The initial `Item` objects as a group.
  - `target`: The final `Item` objects as a group.
  - `trs_keywords`: Keyword arguments for the `Transform` animations.
  - `move_keywords`: Keyword arguments for the `Move` animations.
  - `joint_keywords`: Keyword arguments for the `Succession` animations.

- **Behavior**:
  - It checks that both `source` and `target` are instances of `Group`.
  - It determines whether to perform a move or a direct transform based on the width of the bounding boxes of `source` and `target`.
  - It creates `Succession` animations that combine a move and a transform for each pair of items in `source` and `target`.
  - It returns an `AnimGroup` with the created animations.

### `from_segments` Class Method
This method handles transformations between segments of items.

- **Parameters**:
  - `source`: The initial `Item` object.
  - `source_segments`: An iterable of segment ranges or a single range for the source item.
  - `target`: The final `Item` object or `EllipsisType` to indicate that the target segments should match the source segments.
  - `target_segments`: An iterable of segment ranges or a single range for the target item, or `EllipsisType` to match the source segments.

- **Behavior**:
  - It processes the `source_segments` and `target_segments` to ensure they are in the correct format (lists of tuples representing segment ranges).
  - It creates sub-items for the specified segments in `source` and `target`.
  - It returns a new instance of `TransformInParts` initialized with the segmented `source` and `target` items.

### `matching_patterns` Class Method
This method creates animations based on matching patterns between `source` and `target`.

- **Parameters**:
  - `source`: A tuple containing the initial `MyTypst` object and a list of patterns.
  - `target`: A tuple containing the final `MyTypst` object and a list of patterns.
  - `gapless`: A boolean indicating whether the transformations should be gapless.
  - `from_moved`: A boolean indicating whether the transformations should start from a moved state.
  - `**kwargs`: Additional keyword arguments.

- **Behavior**:
  - It processes the patterns to ensure they are in the correct format.
  - It retrieves the parts of `source` and `target` that match the patterns.
  - Depending on the value of `from_moved`, it either creates `TransformInParts` directly or calls the `from_moved` method.

### `from_transposed` Class Method
This method creates animations based on transposed patterns between `source` and `target`.

- **Parameters**:
  - `source`: The initial `MyTypst` object.
  - `target`: The final `MyTypst` object or `EllipsisType` to indicate that the target patterns should match the source patterns.
  - `source_patterns`: An iterable of patterns for the source.
  - `target_patterns`: An iterable of patterns for the target or `EllipsisType` to match the source patterns.
  - `source_indices`: An iterable of indices for the source patterns or `EllipsisType` to indicate cyclic behavior.
  - `target_indices`: An iterable of indices for the target patterns or `EllipsisType` to indicate cyclic behavior.
  - `**anim_kwargs`: Additional keyword arguments for animations.

- **Behavior**:
  - It processes the input patterns and indices to ensure they are in the correct format.
  - It calls the `matching_patterns` method with the processed patterns and indices.

### Summary
The `TransformInParts` class provides a flexible way to create animations that involve transforming items from a source state to a target state, potentially including intermediate steps like moving. The class supports various types of inputs, including patterns and segment ranges, and allows for customization through keyword arguments.
