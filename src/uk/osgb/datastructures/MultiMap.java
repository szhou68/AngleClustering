package uk.osgb.datastructures;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

public interface MultiMap<K, V> {
	boolean containsKey(K key);
	boolean containsValue(V value);
	Set<K> keySet();
	Collection get(K key);
	boolean put(K key, V value);
	boolean putAll(K key, Collection values);
	Collection remove(K key);
	boolean remove(K key, V value);
	String toString();
	void clear();
	int size();
	Set<Map.Entry<K, V>> entrySet();
}
