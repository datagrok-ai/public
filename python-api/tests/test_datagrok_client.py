from datetime import datetime
import unittest
from unittest.mock import patch, Mock
from datagrok_api import DatagrokClient, Group

class TestDatagrokClient(unittest.TestCase):

    @patch.object(DatagrokClient, '_request')
    def test_get_current_user_group(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {
            'id': 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6',
            'name': 'All users',
            "createdOn": "2025-04-29T09:53:40.502720Z",
            "updatedOn": "2025-04-29T09:53:40.502720Z"
        }
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        group = client.get_current_user_group()
        
        # Assertions to check if the group returned is correct
        self.assertIsInstance(group, Group)
        self.assertEqual(group.id, 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6')
        self.assertEqual(group.name, 'All users')
        self.assertEqual(group.created_on, datetime.fromisoformat("2025-04-29T09:53:40.502720Z"))
        self.assertEqual(group.updated_on, datetime.fromisoformat("2025-04-29T09:53:40.502720Z"))

    @patch.object(DatagrokClient, '_request')
    def test_list_groups(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {'id': 'group1', 'name': 'Test Group 1'},
            {'id': 'group2', 'name': 'Test Group 2'}
        ]
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        groups = client.list_groups()
        
        self.assertEqual(len(groups), 2)
        self.assertEqual(groups[0].id, 'group1')
        self.assertEqual(groups[1].id, 'group2')

    @patch.object(DatagrokClient, '_request')
    def test_lookup_groups(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {'id': 'group1', 'name': 'Test Group'}
        ]
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        groups = client.lookup_groups("Test Group")
        
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0].id, 'group1')
        self.assertEqual(groups[0].name, 'Test Group')

    @patch.object(DatagrokClient, '_request')
    def test_get_group(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = {'id': 'group1', 'name': 'Test Group'}
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        group = client.get_group("group1")
        
        self.assertEqual(group.id, 'group1')
        self.assertEqual(group.name, 'Test Group')

    @patch.object(DatagrokClient, '_request')
    def test_get_group_members(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {'id': 'user1', 'name': 'User One'},
            {'id': 'user2', 'name': 'User Two'}
        ]
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        members = client.get_group_members("group1")
        
        self.assertEqual(len(members), 2)
        self.assertEqual(members[0].id, 'user1')
        self.assertEqual(members[1].id, 'user2')

    @patch.object(DatagrokClient, '_request')
    def test_get_group_memberships(self, mock_request):
        mock_response = Mock()
        mock_response.json.return_value = [
            {'id': 'group1', 'name': 'Test Group'}
        ]
        mock_request.return_value = mock_response
        
        client = DatagrokClient("Bearer token", "http://localhost:8082")
        
        memberships = client.get_group_memberships("group1")
        
        self.assertEqual(len(memberships), 1)
        self.assertEqual(memberships[0].id, 'group1')

if __name__ == '__main__':
    unittest.main()
