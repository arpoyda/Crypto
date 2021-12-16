class Cache:
    def __init__(self, N):
        self.cache = [(-1, -1)] * N

    def __getitem__(self, key):
        index, layer = key
        layer -= 1  # TODO: mb make 1st layer to be 0th
        stored_index = self.cache[layer][0]
        if index == stored_index:
            return self.cache[layer][1]
        return None

    def __setitem__(self, key, value):
        index, layer = key
        layer -= 1
        stored_index = self.cache[layer][0]
        if index != stored_index:  # deny excess assignment
            self.cache[layer] = (index, value)

    def get_auth_path(self):
        return [c[1] for c in self.cache]